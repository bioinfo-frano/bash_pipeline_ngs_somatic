#!/bin/bash

set -euo pipefail

echo "Starting LearnReadOrientationModel (strand artifact modeling)..."

# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"

VARIANT_DIR="$PROJECT_ROOT/data/$SAMPLE/variants"
LOG_DIR="$PROJECT_ROOT/logs/"

JAVA_MEM="-Xmx4g"

F1R2_TAR="$VARIANT_DIR/${SAMPLE}.f1r2.tar.gz"
ORIENTATION_MODEL="$VARIANT_DIR/${SAMPLE}.read-orientation-model.tar.gz"

# ============================================================
# Sanity checks
# ============================================================

echo "Running sanity checks..."

for file in "$F1R2_TAR"
do
  [[ -f "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done

mkdir -p "$LOG_DIR"

echo "All required files found."

# ============================================================
# LearnReadOrientationModel
# ============================================================

echo "Learning read orientation model..."

gatk --java-options "$JAVA_MEM" LearnReadOrientationModel \
  -I "$F1R2_TAR" \
  -O "$ORIENTATION_MODEL" \
  2> "$LOG_DIR/learn_read_orientation_model.log"

echo "LearnReadOrientationModel completed successfully."
echo "Orientation model written to:"
echo "  $ORIENTATION_MODEL"
