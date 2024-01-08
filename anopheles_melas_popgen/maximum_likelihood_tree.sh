# making a maximum likelihood tree containing my melas samples to look at where they fall amongst other species
# include my melas and my anopheles gambiae sensu stricto

# 1. Split the genomic VCF to be mitochondrial only, and then convert the VCF to fasta using --whole genome option
# mitochondrial position on AgamP4_Mt vectorbase 1:15363 

# Create a list of sample names
ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > samples.txt

# Loop through the sample names in samples.txt
#while read -r sample; do
#    # Extract only mitochondrial positions and save as a mitochondrial-only VCF
#    bcftools view -r Mt "${sample}.g.vcf.gz" -O z -o "${sample}_mito_only.g.vcf.gz"
#    tabix -p vcf "${sample}_mito_only.g.vcf.gz"
#
#    # Convert to fasta file
#    python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py \
#        --vcf "${sample}_mito_only.g.vcf.gz" \
#        --ref "/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Mt_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa" \
#        --threads 10 \
#        --snps-no-filt >> "${sample}.log" 2>&1
#done < samples.txt

# I needed to use a multi-sample VCF with the --whole-genome flag

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta_noiupac.py \
        --vcf "mito_only_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz" \
        --ref "/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Mt_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa" \
        --threads 10 \
        --whole-genome

# Take mito fa files for melas from Bijagos and melas from other countries, align this in aliview
# Then add this to holly's aligned fasta file for the tree and realign
# Trim where needed in aliview
# Make tree
raxml-ng --all --msa Sophie_anopheles_mito_2023_aligned.fa --model GTR --prefix Anopheles_mito_1023 --seed 736295 --bs-metric tbe --tree rand{1} --bs-trees 1000

# View in iTOL