#!/bin/bash

set -euo pipefail

TARGET="AGCAGGATAGTCAGTCCTACGTATGGCGGCGCGGTCTATTTGGGGATTGAAC" # this is just 50bp from the centre of the insertion region

OUT="an_coluzzii_specific_region_reads_seqkitfuzzy_50bpregion.tsv"
> "$OUT"

for bam in *.bam; do
  sample=$(basename "$bam" .bam)

  # samtools fasta: filter out secondary (0x100) and supplementary (0x800) with -F 0x900
  # send samtools/seqkit messages to stderr (or silence them with 2>/dev/null if you prefer)
  count=$(
    samtools fasta -F 0x900 "$bam" 2>/dev/null \
    | seqkit grep -s -i -m 2 -p "$TARGET" -n 2>/dev/null \
    | sort -u \
    | wc -l
  )

  echo -e "${sample}\t${count}" >> "$OUT"
done
