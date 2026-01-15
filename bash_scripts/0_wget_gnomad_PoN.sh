#!/bin/bash

set -euo pipefail

echo "Starting downloading: af-only-gnomad.hg38.vcf.gz"

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz

echo "Starting downloading: af-only-gnomad.hg38.vcf.gz.tbi"

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

echo "Starting downloading: 1000g_pon.hg38.vcf.gz"

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz

echo "Starting downloading: 1000g_pon.hg38.vcf.gz.tbi"

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
