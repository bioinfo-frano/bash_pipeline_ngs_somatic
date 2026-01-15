#!/bin/bash

set -euo pipefail

NUMBER_EXPECTED=9

echo "Starting downloading bundle: Homo_sapiens_assembly38"


wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict

wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac
wget -nc https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa


echo "There should be 9 files downloaded"

# Run command in $( )
NUMBER_FILES=$(ls Homo_sapiens_assembly38.* | wc -l)

echo "Number of 'Homo_sapiens_assembly38' files downloaded: $NUMBER_FILES"

# Compare NUMBER_FILES to NUMBER_EXPECTED
if [[ "${NUMBER_FILES}" -eq "${NUMBER_EXPECTED}" ]]; then
  echo "✅ Success: All files were downloaded."
else
  echo "❌ Error: Expected $NUMBER_EXPECTED files but found "$NUMBER_FILES"."
  exit 1
fi
