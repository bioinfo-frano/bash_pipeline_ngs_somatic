#!/bin/bash
set -euo pipefail

echo "Starting GetPileupSummaries..."

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"

REF_FASTA="$PROJECT_ROOT/reference/GRCh38/fasta/Homo_sapiens_assembly38.fasta"
GNOMAD="$PROJECT_ROOT/reference/GRCh38/somatic_resources/af-only-gnomad.hg38.vcf.gz"
INTERVALS="$PROJECT_ROOT/reference/GRCh38/intervals/crc_panel_7genes_sorted.hg38.bed"

BAM="$PROJECT_ROOT/data/$SAMPLE/aligned/${SAMPLE}.sorted.markdup.md.bam"
VARIANT_DIR="$PROJECT_ROOT/data/$SAMPLE/variants"
LOG_DIR="$PROJECT_ROOT/logs"
JAVA_MEM="-Xmx4g"

mkdir -p "$LOG_DIR"

# ============================================================
# GetPileupSummaries
# ============================================================

gatk --java-options "$JAVA_MEM" GetPileupSummaries \
  -R "$REF_FASTA" \
  -I "$BAM" \
  -V "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$VARIANT_DIR/${SAMPLE}.pileups.table" \
  2> "$LOG_DIR/get_pileup_summaries.log"

echo "GetPileupSummaries completed."
