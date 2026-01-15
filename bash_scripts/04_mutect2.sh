#!/bin/bash
set -euo pipefail

echo "Starting somatic variant calling with Mutect2..."

# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"

REF_FASTA="$PROJECT_ROOT/reference/GRCh38/fasta/Homo_sapiens_assembly38.fasta"
REF_FASTA_DICT="$PROJECT_ROOT/reference/GRCh38/fasta/Homo_sapiens_assembly38.dict"
INTERVALS="$PROJECT_ROOT/reference/GRCh38/intervals/crc_panel_7genes_sorted.hg38.bed"

SOMATIC_RESOURCES="$PROJECT_ROOT/reference/GRCh38/somatic_resources"
PON="$SOMATIC_RESOURCES/1000g_pon.hg38.vcf.gz"
GNOMAD="$SOMATIC_RESOURCES/af-only-gnomad.hg38.vcf.gz"

BAM_DIR="$PROJECT_ROOT/data/$SAMPLE/aligned"
INPUT_MD_BAM="$BAM_DIR/${SAMPLE}.sorted.markdup.md.bam"
TUMOR_SM="DMBEL-EIDR-071"                                                           # biological (tumor) sample name

OUTDIR="$PROJECT_ROOT/data/$SAMPLE/variants"
LOG_DIR="$PROJECT_ROOT/logs/"

JAVA_MEM="-Xmx6g"

mkdir -p "$OUTDIR" "$LOG_DIR"

# ============================================================
# Sanity checks (fail early, fail clearly)
# ============================================================

echo "Running sanity checks..."

for file in \
  "$REF_FASTA" \
  "$REF_FASTA_DICT" \
  "$REF_FASTA.fai" \
  "$INPUT_MD_BAM" \
  "$INPUT_MD_BAM.bai" \
  "$INTERVALS" \
  "$PON" \
  "$PON.tbi" \
  "$GNOMAD" \
  "$GNOMAD.tbi"
do
  [[ -f "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done

echo "All required input files found."

# ============================================================
# Step 1: Mutect2 (tumor-only panel of 7 genes)
# ============================================================

echo "Running Mutect2..."

gatk --java-options "$JAVA_MEM" Mutect2 \
  -R "$REF_FASTA" \
  -I "$INPUT_MD_BAM" \
  --tumor-sample "$TUMOR_SM" \
  --panel-of-normals "$PON" \
  --germline-resource "$GNOMAD" \
  -L "$INTERVALS" \
  --af-of-alleles-not-in-resource 0.0000025 \
  --f1r2-tar-gz "$OUTDIR/${SAMPLE}.f1r2.tar.gz" \
  -O "$OUTDIR/${SAMPLE}.unfiltered.vcf.gz" \
  > "$LOG_DIR/mutect2.stdout.log" \
  2> "$LOG_DIR/mutect2.stderr.log"


echo "Mutect2 completed successfully."
echo "Unfiltered VCF written to:"
echo "  $OUTDIR/${SAMPLE}.unfiltered.vcf.gz"
