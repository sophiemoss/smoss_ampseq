# ----------------------------
# The Plan
# Extract read's from each sample's BAM file that map to one of the 4 species amplicons of interest (COX1, IGS, SINE200, MtND4)
# Separate them by sample and amplicon
# Convert them to FASTA
# Blast them to confirm species identity
# Dependencies: bcftools, samtools, seqtk, grep, bash 
# Input files needed: BAM files, reference fasta and species_markers.bed file 
# ----------------------------

# Prepare species markers of interest bed file
# grep -E 'COX1|IGS|SINE200|MtND4' AgamP4_chr.bed > species_markers.bed


#!/bin/bash

# ---- CONFIG ----
BED="/mnt/storage11/sophie/env_dna/species_markers.bed"
REF="/mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa"
THREADS=4

# ---- Process each sample BAM ----
# Index reference if not already
if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference..."
    samtools faidx "$REF"
fi

for BAM in *.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "Processing $SAMPLE..."

    # Ensure BAM is indexed
    [ ! -f "${BAM}.bai" ] && samtools index "$BAM"

    # Loop over each marker region
    while read -r CHR START END NAME; do
        REGION="${CHR}:${START}-${END}"
        echo "  Extracting region $NAME ($REGION)..."

        # 1. Extract region to temp BAM
        samtools view -b "$BAM" "$REGION" > "${SAMPLE}_${NAME}.bam"
        samtools index "${SAMPLE}_${NAME}.bam"

        # 2. Call variants
        bcftools mpileup -f "$REF" "${SAMPLE}_${NAME}.bam" | \
        bcftools call -mv -Oz -o "${SAMPLE}_${NAME}.vcf.gz"

        # 3. Index VCF
        bcftools index "${SAMPLE}_${NAME}.vcf.gz"

        # 4. Generate consensus
        bcftools consensus -f "$REF" "${SAMPLE}_${NAME}.vcf.gz" > "${SAMPLE}_${NAME}_consensus.fasta"

        echo "    → ${SAMPLE}_${NAME}_consensus.fasta created."

    done < "$BED"

    echo "Finished $SAMPLE"
done