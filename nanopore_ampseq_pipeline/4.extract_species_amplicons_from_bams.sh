#!/bin/bash

# ----------------------------
# REQUIRED PACKAGES
# ----------------------------
BED="/mnt/storage11/sophie/env_dna/species_markers.bed"
REF="/mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa"
THREADS=4
OUTDIR="species_markers_merged_fasta"
mkdir -p "$OUTDIR"

# ----------------------------
# Index reference if not already
# ----------------------------
if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference..."
    samtools faidx "$REF"
fi

# ----------------------------
# Process each BAM, extract just the sample name
# ----------------------------
for BAM in *.bam; do
    SAMPLE=$(basename "$BAM" .bam | cut -d_ -f1)
    echo "🧬 Processing $SAMPLE..."

    # Index BAM if needed
    [ ! -f "${BAM}.bai" ] && samtools index "$BAM"

    # Output file for merged consensus
    MERGED_FASTA="${OUTDIR}/${SAMPLE}_species_markers_merged.fasta"
    rm -f "$MERGED_FASTA"

    # Loop through amplicons from BED
    while read -r CHR START END NAME; do
        REGION="${CHR}:${START}-${END}"
        echo "  ↪ Extracting $NAME from $REGION..."

        # Temp files
        TEMP_BAM="${SAMPLE}_${NAME}.bam"
        VCF_GZ="${SAMPLE}_${NAME}.vcf.gz"

        # Extract BAM region
        samtools view -b "$BAM" "$REGION" > "$TEMP_BAM"

        # Skip if BAM is empty
        if [ ! -s "$TEMP_BAM" ]; then
            echo "    ⚠️  No reads in $TEMP_BAM, skipping $NAME."
            rm -f "$TEMP_BAM"
            continue
        fi

        samtools index "$TEMP_BAM"

        # Variant calling
        bcftools mpileup -f "$REF" "$TEMP_BAM" | \
        bcftools call -mv -Oz -o "$VCF_GZ"

        bcftools index -f "$VCF_GZ"

        # Extract consensus sequence for the specific amplicon region
        bcftools consensus -f "$REF" "$VCF_GZ" > tmp_consensus.fa
        samtools faidx tmp_consensus.fa "$REGION" > tmp_region.fa
        SEQ=$(grep -v "^>" tmp_region.fa | tr -d '\n')

        if [ -z "$SEQ" ]; then
            echo "    ⚠️  Empty sequence for $NAME — skipping."
        else
            echo ">${SAMPLE}_${NAME}" >> "$MERGED_FASTA"
            echo "$SEQ" >> "$MERGED_FASTA"
            echo "    ✔ $NAME added to ${MERGED_FASTA}"
        fi

        # Cleanup
        rm -f "$TEMP_BAM" "$TEMP_BAM.bai" "$VCF_GZ" "$VCF_GZ.csi"

    done < "$BED"

    echo "✅ Finished $SAMPLE → ${MERGED_FASTA}"
done

echo "🎯 All merged FASTAs are saved in: $OUTDIR/"
