# Using the dorado pipeline for basecalling and aligning env data to anopheles gambiae

dorado_pipe \
dorado --sample_sheet /mnt/storageG1/data/experiments/Exp149_eDNA_amplicon_Ghana/exp149/20250304_1604_MN43600_FBA90131_4542a511/exp149.csv \
--barcoding_kit SQK-NBD114-96 \
--ID Exp149_eDNA_amplicon_Ghana \
--data_directory /mnt/storageG1/data/experiments/Exp149_eDNA_amplicon_Ghana/exp149/  \
--output_directory /mnt/storageG1/data/experiments/Exp149_eDNA_amplicon_Ghana/exp149/rerun_nanopore_pipeline_out_exp149/ \
--basecalling_model dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
--ref_seq /mnt/storageG1/data/experiments/Exp149_eDNA_amplicon_Ghana/reference_genome/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
--kraken true \
--threads 20 \
--trimmer none


#2 run the extra demultiplexing (I did this in G1)

python demux_nanopore_plates_sophie.py \
-p plate_layout.csv \
-b internal_barcodes.csv \
-f /mnt/storageG1/data/experiments/Exp149_eDNA_amplicon_Ghana/exp149/nanopore_pipeline_out_exp149/fastq/fastq \
-o samples.csv \
-m 0 \
-t 8

#3 Move the fastq files that have been double demultiplexed to s11 and then continue with amplicon script.

ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv

#4
for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done

#5 Run the amplicons script


# nanopore_amplicon_script.py [-h] --index-file INDEX_FILE --ref REF --gff GFF --bed BED
#                                             [--min-base-qual MIN_BASE_QUAL]
#                                             [--version]
#
#Nanopore amplicon sequencing analysis script
#
#optional arguments:
#  -h, --help            show this help message and exit
#  --index-file INDEX_FILE
#                        sample_file.csv with the "sample" column for sample IDs (default: None)
#  --ref REF             Reference fasta (default: None)
#  --gff GFF             GFF file (default: None)
#  --bed BED             BED file with genes/amplicon locations (default: None)
#  --min-base-qual MIN_BASE_QUAL
#                        Minimum base quality to use by freebayes (default: 30)
#  --version             show program's version number and exit

#3

ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv

#4
for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done

#5

python /mnt/storage11/sophie/gitrepos/amplicon-seq/scripts/sophie_nanopore_amplicon_script.py \
--index-file sample_file.csv \
--ref /mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
--gff /mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.56.gff3 \
--bed /mnt/storage11/sophie/env_dna/AgamP4_chr.bed \
--min-base-qual 30 \
--threads 40 > amplicon_log.txt 2>&1

# UNKN in bed file? will this work to find IGS and SINE200 and ITS1 and ITS2? Should do.

# Checks on correct bed file before running:

# Print any region in the bed file that does not match with an annotation in the gff file:
bedtools intersect -a AgamP4_chr.bed -b Anopheles_gambiae.AgamP4.56.gff3 -wa | sort | uniq > matched.bed
comm -23 <(cut -f1-4 AgamP4_chr.bed | sort) matched.bed
# this prints all regions that match between the bed file and the gff file to matched.bed and then prints any 
# that do not match in the terminal for you to see.

# Check that regions in the bedfile match the reference fastq file
awk 'NR==FNR {chrlen[$1]=$2; next} $1 in chrlen && $3 <= chrlen[$1] {next} {print "Invalid:", $0}' Anopheles_gambiae.AgamP4.dna.toplevel.fa.fai AgamP4_chr.bed
# If nothing prints, you're good.

#6 Now run the amplicon script

## running with:
python /mnt/storage11/sophie/gitrepos/amplicon-seq/scripts/sophie_nanopore_amplicon_script_minimap2_nosamclip.py \
--index-file sample_file.csv --ref /mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
--gff /mnt/storage11/sophie/env_dna/Anopheles_gambiae.AgamP4.56.gff3 --bed /mnt/storage11/sophie/env_dna/AgamP4_chr.bed \
--min-base-qual 30 --threads 20 > amplicon_log.txt 2>&1

#7 Create amplicon x sample coverage matrix by combining coverage files (the sample_coverage_mean.txt files) to see how good they are across the sample set




