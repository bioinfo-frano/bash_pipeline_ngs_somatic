#!/bin/bash
set -euo pipefail

echo "Starting CalculateContamination..."

# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"

VARIANT_DIR="$PROJECT_ROOT/data/$SAMPLE/variants"
LOG_DIR="$PROJECT_ROOT/logs"
JAVA_MEM="-Xmx2g"

PILEUP_TABLE="$VARIANT_DIR/${SAMPLE}.pileups.table"
CONTAM_TABLE="$VARIANT_DIR/${SAMPLE}.contamination.table"

mkdir -p "$LOG_DIR"

# ============================================================
# Sanity checks
# ============================================================

for file in "$PILEUP_TABLE"
do
  [[ -f "$file" ]] || { echo "ERROR: Missing $file"; exit 1; }
done

# ============================================================
# CalculateContamination
# ============================================================

gatk --java-options "$JAVA_MEM" CalculateContamination \
  -I "$PILEUP_TABLE" \
  -O "$CONTAM_TABLE" \
  2> "$LOG_DIR/calculate_contamination.log"

echo "CalculateContamination completed."
echo "Contamination table written to:"
echo "  $CONTAM_TABLE"
