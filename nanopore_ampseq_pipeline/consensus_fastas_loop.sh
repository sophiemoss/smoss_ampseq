#!/bin/bash
set -euo pipefail

# Set reference genome
REFERENCE="Anopheles_gambiae.AgamP4.dna.toplevel.fa"

# Index reference if needed
if [ ! -f "$REFERENCE.fai" ]; then
    echo "Indexing reference FASTA..."
    samtools faidx "$REFERENCE"
fi

# Loop through each BAM file
for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    echo "🧬 Processing sample: $sample"

    # Step 1: Variant calling with GATK (GVCF mode)
    gatk HaplotypeCaller \
        -I "$bam" \
        -R "$REFERENCE" \
        -O "$sample.vcf.gz" \
        -ERC GVCF

    # Step 2: Create mask BED of low-depth sites (DP < 10)
    bcftools view -i 'FMT/DP<10' "$sample.vcf.gz" \
    | bcftools convert --gvcf2vcf -f "$REFERENCE" \
    | bcftools query -f '%CHROM\t%POS\t%POS\n' > "$sample.mask.bed"

    # Step 3: Extract SNPs only (no indels)
    bcftools view -V indels -c 1 -a "$sample.vcf.gz" -Oz -o "$sample.snps.vcf.gz"
    tabix -p vcf "$sample.snps.vcf.gz"

    # Step 4: Generate consensus FASTA with mask, replace * with N
    bcftools consensus -f "$REFERENCE" -m "$sample.mask.bed" "$sample.snps.vcf.gz" \
    | tr '*' 'N' > "$sample.consensus.fasta"

    samtools faidx "$sample.consensus.fasta"

    # Step 5: Extract specific regions and rename headers to include sample and amplicon name
    samtools faidx "$sample.consensus.fasta" Mt:1425-2960 \
        | sed "s/^>.*/>${sample}_COX1/" > "${sample}.COX1.fa"

    samtools faidx "$sample.consensus.fasta" UNKN:36384025-36384413 \
        | sed "s/^>.*/>${sample}_IGS/" > "${sample}.IGS.fa"

    samtools faidx "$sample.consensus.fasta" X:22951332-22951809 \
        | sed "s/^>.*/>${sample}_SINE200/" > "${sample}.SINE200.fa"

    samtools faidx "$sample.consensus.fasta" UNKN:31004851-31005267 \
        | sed "s/^>.*/>${sample}_ITS1/" > "${sample}.ITS1.fa"

    samtools faidx "$sample.consensus.fasta" UNKN:35962838-35963328 \
        | sed "s/^>.*/>${sample}_ITS2/" > "${sample}.ITS2.fa"

    samtools faidx "$sample.consensus.fasta" Mt:8309-8824 \
        | sed "s/^>.*/>${sample}_MtND4/" > "${sample}.MtND4.fa"

    echo "✅ Finished processing $sample"
done

# 🧬 Step 6: Concatenate all per-region FASTAs into combined files
echo "📦 Concatenating per-region FASTAs..."

cat *.COX1.fa    > combined.COX1.fasta
cat *.IGS.fa     > combined.IGS.fasta
cat *.SINE200.fa > combined.SINE200.fasta
cat *.ITS1.fa    > combined.ITS1.fasta
cat *.ITS2.fa    > combined.ITS2.fasta
cat *.MtND4.fa   > combined.MtND4.fasta

echo "✅ Combined FASTA files created:"
echo " - combined.COX1.fasta"
echo " - combined.IGS.fasta"
echo " - combined.SINE200.fasta"
echo " - combined.ITS1.fasta"
echo " - combined.ITS2.fasta"
echo " - combined.MtND4.fasta"
