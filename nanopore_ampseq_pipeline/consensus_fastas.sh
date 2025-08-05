# Create individual vcf files for each sample, variant call with gatk (MAKE SURE USING gatk=4.1.4.1) ##### VERY IMPORTANT #####

gatk HaplotypeCaller -I ADAR210.bam -R Anopheles_gambiae.AgamP4.dna.toplevel.fa -O ADAR210.vcf.gz -ERC GVCF

# Now create a mask over all positions with DP less than 10 (for diploid)

bcftools view -i 'FMT/DP<10' ADAR210.vcf.gz | bcftools convert --gvcf2vcf -f Anopheles_gambiae.AgamP4.dna.toplevel.fa| bcftools query -f '%CHROM\t%POS\t%POS\n' > ADAR210.mask.bed

# Create a VCF containing only SNPs (not indels)

bcftools view -V indels -c 1 -a ADAR210.vcf.gz -Oz -o ADAR210.snps.vcf.gz
tabix -p vcf ADAR210.snps.vcf.gz

# Use the filtered VCF and the mas bed file to create the consensus sequence

bcftools consensus -f Anopheles_gambiae.AgamP4.dna.toplevel.fa -m ADAR210.mask.bed ADAR210.snps.vcf.gz | tr '*' 'N' > ADAR210.consensus.fasta
samtools faidx ADAR210.consensus.fasta


# Split the consensus fasta to different regions
samtools faidx ADAR210.consensus.fasta Mt:1425-2960 | sed 's/^>.*/>ADAR210_COX1/' > ADAR210.COX1.fa
samtools faidx ADAR210.consensus.fasta UNKN:36384025-36384413 | sed 's/^>.*/>ADAR210_IGS/' > ADAR210.IGS.fa
samtools faidx ADAR210.consensus.fasta X:22951332-22951809 | sed 's/^>.*/>ADAR210_SINE200/' > ADAR210.SINE200.fa
samtools faidx ADAR210.consensus.fasta UNKN:31004851-31005267 | sed 's/^>.*/>ADAR210_ITS1/' > ADAR210.ITS1.fa
samtools faidx ADAR210.consensus.fasta UNKN:35962838-35963328 | sed 's/^>.*/>ADAR210_ITS2/' > ADAR210.ITS2.fa
samtools faidx ADAR210.consensus.fasta Mt:8309-8824 | sed 's/^>.*/>ADAR210_MtND4/' > ADAR210.MtND4.fa



