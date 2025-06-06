# Using the dorado pipeline for basecalling (multiple species so no alignment to ref genomes)

dorado basecaller \
--models-directory /mnt/storageG1/data/dorado/models/dna/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
--kit-name SQK-NBD114-96 dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/ > /mnt/storageG1/data/experiments/Exp178_AEDES_AP_AE_AMP/calls/calls.bam

#dorado demux step to produce the individual fastq per barcode 
dorado demux --output-dir ./fastq/ --emit-fastq --no-classify ./calls/calls.bam

#2 Move the fastq files that have been demultiplexed to s11 and then continue with amplicon script.

# zip all fastq files

for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done

# create sample_file.csv

ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv

# 3 Run checks before amplicon pipeline

# Print any region in the bed file that does not match with an annotation in the gff file:
bedtools intersect -a aedes_minion_renamed_chrs.bed -b GCF_002204515.2_AaegL5.0_genomic.gff -wa | sort | uniq > matched.bed
comm -23 <(cut -f1-4 aedes_minion_renamed_chrs.bed | sort) matched.bed


###_______________________________________________________________________###
###__________ 4 Run amplicon pipeline for Aedes aegypti samples___________###
###_______________________________________________________________________###


# Make sure there is no samclip as this gets rid of all amplicon reads
# This script takes a SINGLE FASTQ per sample. Use minimap2 and allow freebayes to chunk genome.

python /mnt/storage11/sophie/gitrepos/anopheles_ampseq/nanopore_ampseq_pipeline/sophie_nanopore_amplicon_script_minimap2.py \
--index-file sample_file.csv \
--ref GCF_002204515.2_AaegL5.0_genomic.fna \
--gff GCF_002204515.2_AaegL5.0_genomic.gff \
--bed aedes_minion_renamed_chrs.bed \
--min-base-qual 20 \
--threads 40 > amplicon_log.txt 2>&1

## Now create coverage matrix to show coverage of amplicons

python /mnt/storage11/sophie/gitrepos/anopheles_ampseq/nanopore_ampseq_pipeline/2.create_ampiconxsample_coverage_matrix.py

## Now compute SNP allele frequencies taking into account missingness and filtering for 20 so it's more stringent

python /mnt/storage11/sophie/gitrepos/anopheles_ampseq/nanopore_ampseq_pipeline/3.compute_SNP_allele_frequencies.py



###__________________________________________________________________________###
###__________ 4 Run amplicon pipeline for Aedes albopictus samples___________###
###__________________________________________________________________________###

python script.py \
--index-file sample_file.csv \ 
--ref GCF_035046485.1_AalbF5_genomic.fna \
--gff GCF_035046485.1_AalbF5_genomic.gff \
--bed aedes_minion.bed \
--min-base-qual 30 \
--threads 40 > amplicon_log.txt 2>&1

# Albopictus chromosomes
NC_085136.1
NC_085137.1
NC_085138.1
NC_006817.1

# 7 Create amplicon x sample coverage matrix by combining coverage files (the sample_coverage_mean.txt files) to see how good they are across the sample set




