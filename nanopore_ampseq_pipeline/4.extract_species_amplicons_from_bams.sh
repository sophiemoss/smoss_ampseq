# ----------------------------
# The Plan
# Extract read's from each sample's BAM file that map to one of the 4 species amplicons of interest (COX1, IGS, SINE200, MtND4)
# Separate them by sample and amplicon
# Convert them to FASTA
# Blast them to confirm species identity
# ----------------------------

# Prepare regions of interest bed file
# grep -E 'COX1|IGS|SINE200|MtND4' AgamP4_chr.bed > species_markers.bed


#!/bin/bash

# ---- CONFIG ----
BED="/mnt/storage11/sophie/env_dna/species_markers.bed"
THREADS=4

# ---- Process each sample BAM ----
for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    echo "Processing $sample..."

    # Ensure BAM is indexed
    if [ ! -f "${bam}.bai" ]; then
        echo "Indexing $bam..."
        samtools index "$bam"
    fi

    # Loop through each marker amplicon region
    while read -r chrom start end name; do
        region="${chrom}:${start}-${end}"
        echo "  Extracting $name from $region..."

        # Extract reads from region, convert to FASTA
        samtools view -b "$bam" "$region" | \
            samtools fastq - | \
            seqtk seq -A - > "${sample}_${name}.fasta"

    done < "$BED"

    echo "Done with $sample"
done
