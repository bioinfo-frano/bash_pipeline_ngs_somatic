#!/bin/bash
set -euo pipefail

# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"

VARIANT_DIR="$PROJECT_ROOT/data/$SAMPLE/variants"
LOG_DIR="$PROJECT_ROOT/logs"
mkdir -p "$LOG_DIR"

FILTERED_VCF="$VARIANT_DIR/${SAMPLE}.filtered.vcf.gz"
POSTFILTER_VCF="$VARIANT_DIR/${SAMPLE}.postfiltered.vcf.gz"
SUMMARY_TXT="$VARIANT_DIR/${SAMPLE}.postfilter_summary.txt"
LOG_FILE="$LOG_DIR/${SAMPLE}.postfilter.log"

# ============================================================
# Thresholds (amplicon tumor-only)
# ============================================================

MIN_DP=200          # total depth
MIN_AD_ALT=10       # ALT read count
MIN_VAF=0.02        # 2%

# ============================================================
# Logging
# ============================================================

exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date)] Starting post-filtering"
echo "Sample: $SAMPLE"

# ============================================================
# Sanity checks
# ============================================================

for file in "$FILTERED_VCF" "$FILTERED_VCF.tbi"
do
  [[ -f "$file" ]] || { echo "ERROR: Missing $file"; exit 1; }
done

# ============================================================
# Apply hard filters
# ============================================================

echo "Applying post-filter thresholds:"
echo "  DP >= $MIN_DP"
echo "  ALT reads (AD[1]) >= $MIN_AD_ALT"
echo "  VAF >= $MIN_VAF"

bcftools view -f PASS "$FILTERED_VCF" | \
bcftools filter \
  -i "FORMAT/DP >= ${MIN_DP} && FORMAT/AD[0:1] >= ${MIN_AD_ALT} && FORMAT/AF >= ${MIN_VAF}" \
  -Oz -o "$POSTFILTER_VCF"

# ============================================================
# Index the post-filtered VCF (required for IGV) → .tbi
# ============================================================

echo "Indexing post-filtered VCF"

bcftools index -t "$POSTFILTER_VCF" # Option '-t' → .tbi. Without any option → .csi (default). Both are valid index.

# Sanity check: ensure index was created
if [[ ! -f "${POSTFILTER_VCF}.tbi" ]]; then
  echo "ERROR: Tabix index (.tbi) was not created"
  exit 1
fi

# ============================================================
# Variant counts
# ============================================================

N_VARIANTS=$(bcftools view -H "$POSTFILTER_VCF" | wc -l)

if [[ "$N_VARIANTS" -eq 0 ]]; then
  echo "WARNING: 0 variants passed post-filtering"
else
  echo "Variants retained after post-filtering: $N_VARIANTS"
fi

# ============================================================
# Summary file
# ============================================================

{
  echo "Post-filter summary"
  echo "========================"
  echo "Sample: $SAMPLE"
  echo ""
  echo "Library type: Amplicon (PCR)"
  echo "Sequencing: Tumor-only"
  echo ""
  echo "Thresholds:"
  echo "  DP >= $MIN_DP"
  echo "  ALT reads >= $MIN_AD_ALT"
  echo "  VAF >= $MIN_VAF"
  echo ""
  echo "Variants retained: $N_VARIANTS"
} > "$SUMMARY_TXT"

echo "[$(date)] $SAMPLE post-filtering completed successfully"
