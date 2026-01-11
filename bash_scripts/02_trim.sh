#!/bin/bash

set -euo pipefail

echo "Starting Cutadapt..."

# --- Configuration ---
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRA_ID"

RAW_FASTQ_DIR="$PROJECT_ROOT/data/$SAMPLE/raw_fastq"
TRIM_DIR="$PROJECT_ROOT/data/$SAMPLE/trimmed"
QC_DIR="$PROJECT_ROOT/data/$SAMPLE/qc"
LOG_DIR="$PROJECT_ROOT/logs/"

THREADS=4

# --- Create output directory ---
mkdir -p "$TRIM_DIR"
mkdir -p "$QC_DIR"
mkdir -p "$LOG_DIR"

# --- Run Cutadapt ---
# -u 5 -u -5:  Trim 5 bp from 5′ and 3' ends of R1
# -U 5 -U -5:  Trim 5 bp from 5′ and 3' ends of R2
# -m 30: Discard reads shorter than 30 bp after trimming
# -q 20,20: Trims low-quality bases from 5' and 3′ ends of reads. 20 corresponds to Phred Q20 for 5' and 3' ends, respectively. Trimms only when base quality is < Q20. Applies to both R1 and R2.
# -a A{10}: Trims poly-A tails of ≥10 consecutive A from the 3′ end of R1 (Internal A-rich regions are not affected. Only tails with 10 or more consecutive A’s are trimmed)
# -A A{10}: Trims poly-A tails of ≥10 consecutive A from the 3′ end of R2 (Internal A-rich regions are not affected. Only tails with 10 or more consecutive A’s are trimmed)

# --report=full Generate full report
cutadapt \
  -u 5 -u -5 \
  -U 5 -U -5 \
  -q 20,20 \
  -m 30 \
  -a A{10} \
  -A A{10} \
  -j "$THREADS" \
  --report=full \
  -o "$TRIM_DIR/SRR30536566_R1.trimmed.fastq.gz" \
  -p "$TRIM_DIR/SRR30536566_R2.trimmed.fastq.gz" \
  "$RAW_FASTQ_DIR/SRR30536566_1.fastq.gz" \
  "$RAW_FASTQ_DIR/SRR30536566_2.fastq.gz" \
  > "$LOG_DIR/cutadapt_SRR30536566.log"

echo "Cutadapt completed successfully."

# --- Run FastQC ---
fastqc \
  --threads "$THREADS" \
  --outdir "$QC_DIR" \
  "$TRIM_DIR"/*.fastq.gz

echo "FastQC completed successfully."

# --- Run MultiQC ---
multiqc "$QC_DIR" --outdir "$QC_DIR"

echo "MultiQC completed successfully."
