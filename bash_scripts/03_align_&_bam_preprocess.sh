#!/bin/bash

set -euo pipefail

echo "Starting alignment and BAM preprocessing..."

# =========================
# Configuration
# =========================

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"

REF_DIR="$PROJECT_ROOT/reference/GRCh38/fasta"
REF_FASTA="$REF_DIR/Homo_sapiens_assembly38.fasta"

TRIM_DIR="$PROJECT_ROOT/data/$SAMPLE/trimmed"
ALIGN_DIR="$PROJECT_ROOT/data/$SAMPLE/aligned"
QC_DIR="$PROJECT_ROOT/data/$SAMPLE/qc"
LOG_DIR="$PROJECT_ROOT/logs/"
FINAL_BAM="$ALIGN_DIR/${SAMPLE}.sorted.markdup.md.bam"

THREADS=4

# Read group information (REQUIRED by GATK)
RG_ID="SRR30536566"
RG_SM="DMBEL-EIDR-071"
RG_LB="AMPLICON"
RG_PL="ILLUMINA"
RG_PU="HiSeq4000"

# =========================
# Create directories
# =========================

mkdir -p "$ALIGN_DIR" "$QC_DIR" "$LOG_DIR"

# =========================
# Step 0: Check / build BWA index
# =========================

echo "Checking BWA index..."

if [[ ! -f "${REF_FASTA}.64.bwt" && ! -f "${REF_FASTA}.bwt" ]]; then
  echo "BWA index not found for reference genome."
  echo "Building BWA index for reference genome..."
  bwa index "$REF_FASTA" 2> "$LOG_DIR/bwa_index.log"
  echo "BWA indexing completed."
else
  echo "BWA index found. Skipping indexing of reference genome."
fi

# =========================
# Step 1: Alignment: BWA-MEM with Read Groups → SAM
# =========================

echo "Running BWA-MEM alignment..."

bwa mem \
  -t "$THREADS" \
  -R "@RG\tID:${RG_ID}\tSM:${RG_SM}\tLB:${RG_LB}\tPL:${RG_PL}\tPU:${RG_PU}" \
  "$REF_FASTA" \
  "$TRIM_DIR/SRR30536566_R1.trimmed.fastq.gz" \
  "$TRIM_DIR/SRR30536566_R2.trimmed.fastq.gz" \
  > "$ALIGN_DIR/${SAMPLE}.sam" \
  2> "$LOG_DIR/bwa_mem.log"

echo "Alignment completed."

# =========================
# Step 2: Convert SAM → BAM
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE}.sam" ]]; then
  echo "ERROR: SAM file not created!" >&2
  exit 1
fi

echo "Converting SAM to BAM..."

samtools view -@ "$THREADS" -bS \
  "$ALIGN_DIR/${SAMPLE}.sam" \
  > "$ALIGN_DIR/${SAMPLE}.bam"

# =========================
# Step 3: Sort BAM (coordinate sort)
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE}.bam" ]]; then
  echo "ERROR: BAM file not created!" >&2
  exit 1
fi

echo "Sorting BAM..."

samtools sort -@ "$THREADS" \
  "$ALIGN_DIR/${SAMPLE}.bam" \
  "$ALIGN_DIR/${SAMPLE}.sorted"

rm "$ALIGN_DIR/${SAMPLE}.sam" "$ALIGN_DIR/${SAMPLE}.bam"

echo "Sorting completed."

# =========================
# Step 4: Mark duplicates (AMPLICON-AWARE: TAG DUPLICATES)
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE}.sorted.bam" ]]; then
  echo "ERROR: Sorted BAM not created!" >&2
  exit 1
fi

echo "Marking duplicates..."

picard MarkDuplicates \
  INPUT="$ALIGN_DIR/${SAMPLE}.sorted.bam" \
  OUTPUT="$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" \
  METRICS_FILE="$ALIGN_DIR/${SAMPLE}.markdup.metrics.txt" \
  CREATE_INDEX=false \
  REMOVE_DUPLICATES=false \
  TAG_DUPLICATE_SET_MEMBERS=true \
  VALIDATION_STRINGENCY=SILENT \
  2> "$LOG_DIR/markduplicates.log"

echo "Duplicate marking completed."

# =========================
# Step 5: Add MD and NM tags (GATK robustness) + Index final BAM (REQUIRED for GATK)
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" ]]; then
  echo "ERROR: MarkDuplicates failed!" >&2
  exit 1
fi

echo "Adding MD tags..."

samtools calmd -b \
  "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" \
  "$REF_FASTA" \
  > "$FINAL_BAM"

if [[ ! -s "$FINAL_BAM" ]]; then
  echo "ERROR: Final BAM not created by samtools calmd!" >&2
  exit 1
fi

echo "Indexing final BAM..."

samtools index "$FINAL_BAM"

echo "BAM indexing completed."

# =========================
# Step 6: Alignment statistics
# =========================

echo "Generating alignment statistics..."

samtools flagstat \
  "$FINAL_BAM" \
  > "$LOG_DIR/${SAMPLE}.flagstat.txt"

# =========================
# Step 7: MultiQC (duplicates + alignment metrics)
# =========================

echo "Running MultiQC..."

multiqc \
  "$ALIGN_DIR" \
  "$LOG_DIR" \
  --outdir "$QC_DIR"

# =========================
# Step 8: Cleanup intermediate files (Optional)
# =========================

echo "Cleaning up intermediate files..."

rm -f \
  "$ALIGN_DIR/${SAMPLE}.sorted.bam" \
  "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam"

echo "Cleanup completed."


echo "Alignment and BAM preprocessing completed successfully."
