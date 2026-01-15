#!/bin/bash
set -euo pipefail

# =========================
# Configuration
# =========================
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GTF="$PROJECT_ROOT/reference/GRCh38/intervals/gencode.v38.annotation.gtf"
OUT_BED="$PROJECT_ROOT/reference/GRCh38/intervals/crc_panel_7genes.hg38.bed"
OUT_BED_SORTED="$PROJECT_ROOT/reference/GRCh38/intervals/crc_panel_7genes_sorted.hg38.bed"

echo "Generating BED file for CRC 7-gene panel..."

# =========================
# Step 1: Extract gene intervals from GTF
# =========================

awk '
BEGIN {
  FS="\t"; OFS="\t"
}

# Skip header/comment lines
/^#/ { next }

# Process only gene features
$3=="gene" {
  gene=""
  n=split($9, attrs, ";")
  for (i=1; i<=n; i++) {
    if (attrs[i] ~ /gene_name "/) {
      gsub(/.*gene_name "|"/, "", attrs[i])
      gene=attrs[i]
    }
  }

  # Keep only CRC-relevant genes
  if (gene=="KRAS" || gene=="NRAS" || gene=="BRAF" ||
      gene=="PIK3CA" || gene=="PTEN" || gene=="RRAS" ||
      gene=="MAP2K1") {
    # BED format is 0-based start, 1-based end
    print $1, $4-1, $5, gene
  }
}
' "$GTF" > "$OUT_BED"

echo "Unsorted BED file created:"
ls -lh "$OUT_BED"

# =========================
# Step 2: Sort BED file (recommended for GATK)
# =========================

sort -k1,1 -k2,2n "$OUT_BED" -o "$OUT_BED_SORTED"

echo "Sorted BED file created:"
ls -lh "$OUT_BED_SORTED"

# =========================
# Step 3: Display final BED
# =========================

echo
echo "Final CRC 7-gene BED file contents:"
cat "$OUT_BED_SORTED"
