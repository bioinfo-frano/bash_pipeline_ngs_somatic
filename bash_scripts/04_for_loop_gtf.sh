#!/bin/bash

GTF="gencode.v38.annotation.gtf.gz"
GENES=("KRAS" "NRAS" "BRAF" "PIK3CA" "PTEN" "RRAS" "MAP2K1")

echo "Checking for required genes in $GTF..."

for gene in "${GENES[@]}"
do
  echo -n "Checking $gene... "

  # Capture the first matching line for this gene
  match=$(zless "$GTF" | grep "^chr" | grep -w "$gene" | head -n1)

  # Check if we found a match
  if [ -n "$match" ]
  then
    echo "FOUND"
    echo "  Evidence: $match"
  else
    echo "NOT FOUND"
    echo "ERROR: Missing required gene: $gene"
    exit 1
  fi
  echo ""
done

echo "All genes were found in $GTF"
