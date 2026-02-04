#!/bin/bash

# Detects region shared in both An. coluzzii and An. gambiae, close to the coluzzii SINE200 X6.1 insertion.
# This code scans through BAM files and counts how many sequencing reads contian the specific TARGET DNA sequence.
# For each BAM, the script:
# Extracts all sequencing reads as FASTA sequences
# Searches for reads that match a given DNA target sequence (allowing fuzzy matches which are specified with -m)
# Counts how many unique read names match
# Writes the results into a tsv

set -euo pipefail

TARGET="ACATCCTAAAATAATGAATTAAATGCAGTTCCTTATTATTCTCAGCTAACCGCTATTGTTCAGCGTGAACAATGCATTAAATTCTGAGTGGAGTAAAATGTATC"
OUT="m2_an_colgamb_region_reads_seqkitfuzzy.txt"
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


