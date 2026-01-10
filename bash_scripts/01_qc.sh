#!/bin/bash

set -euo pipefail

echo "Starting FastQC..."

# --- Configuration ---
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRA_ID"

RAW_FASTQ_DIR="$PROJECT_ROOT/data/$SAMPLE/raw_fastq"
QC_DIR="$PROJECT_ROOT/data/$SAMPLE/qc"

THREADS=4

# --- Create output directory ---
mkdir -p "$QC_DIR"

# --- Run FastQC ---
fastqc \
  --threads "$THREADS" \
  --outdir "$QC_DIR" \
  "$RAW_FASTQ_DIR"/*.fastq.gz

echo "FastQC completed successfully."

# --- Run MultiQC ---
multiqc "$QC_DIR" --outdir "$QC_DIR"

echo "MultiQC completed successfully."
