#!/bin/bash
set -euo pipefail

# empty the output files at the start


for bam in *.bam; do
    count=$(samtools view "$bam" | grep -c 'TTCAGTGACTAATTGATTCCCCCCCATAGCAGGATAGTCAGTCCTACGTATGGCGGCGCGGTCTATTTGGGGATTGAACCCATGACAGGCATGTTGTTAAG')
    sample=$(basename "$bam" .bam)

    echo -e "${sample}\t${count}" >> grep_an_coluzzii_specific_region_reads.txt
done
