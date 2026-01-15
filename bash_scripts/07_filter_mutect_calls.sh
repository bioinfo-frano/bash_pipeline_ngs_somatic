#!/bin/bash

set -euo pipefail

echo "Starting FilterMutectCalls (final somatic filtering)..."

# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"

REF_FASTA="$PROJECT_ROOT/reference/GRCh38/fasta/Homo_sapiens_assembly38.fasta"

VARIANT_DIR="$PROJECT_ROOT/data/$SAMPLE/variants"
LOG_DIR="$PROJECT_ROOT/logs"

JAVA_MEM="-Xmx4g"

UNFILTERED_VCF="$VARIANT_DIR/${SAMPLE}.unfiltered.vcf.gz"
ORIENTATION_MODEL="$VARIANT_DIR/${SAMPLE}.read-orientation-model.tar.gz"
FILTERED_VCF="$VARIANT_DIR/${SAMPLE}.filtered.vcf.gz"
CONTAM_TABLE="$VARIANT_DIR/${SAMPLE}.contamination.table"

# ============================================================
# Sanity checks
# ============================================================

echo "Running sanity checks..."

for file in \
  "$REF_FASTA" \
  "$UNFILTERED_VCF" \
  "$UNFILTERED_VCF.tbi" \
  "$ORIENTATION_MODEL" \
  "$CONTAM_TABLE"
do
  [[ -f "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done

mkdir -p "$LOG_DIR"

echo "All required files found."

# ============================================================
# FilterMutectCalls
# ============================================================

echo "Filtering Mutect2 calls..."

gatk --java-options "$JAVA_MEM" FilterMutectCalls \
  -R "$REF_FASTA" \
  -V "$UNFILTERED_VCF" \
  --contamination-table "$CONTAM_TABLE" \
  --orientation-bias-artifact-priors "$ORIENTATION_MODEL" \
  -O "$FILTERED_VCF" \
  2> "$LOG_DIR/filter_mutect_calls.log"

echo "FilterMutectCalls completed successfully."
echo "Final filtered VCF:"
echo "  $FILTERED_VCF"
