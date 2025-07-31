#1 Basecalling the data from the minion, including demux for the ONT barcodes
dorado basecaller \
--models-directory /mnt/storageG1/data/dorado/models/dna/dna_r10.4.1_e8.2_400bps_sup@v5.0.0/ \
--kit-name SQK-NBD114-96 sup /mnt/storageG1/data/experiments/Exp195_kai/Exp195/20250728_1644_MN43600_FBC21767_37bd1382/pod5 > /mnt/storageG1/data/experiments/Exp195_kai/Exp195/nanopore_pipeline_out/calls.bam
 
 
dorado demux --output-dir ./fastq/ --emit-fastq --no-classify ./calls/calls.bam

#2 Run the extra demultiplexing to sort between our internal 380 barcodes
# make the plate_layout.csv and the internal_barcodes.csv
# the samples.csv is an output file

for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done

# Rename fastq files
# Loop over matching files in the current directory
for file in *_barcode*.fastq*; do
  # Extract just the barcode and extension
  newname=$(echo "$file" | sed -E 's/.*_(barcode[0-9]+\.fastq(.gz)?)$/\1/')
  # Rename the file if the pattern matched
  if [[ "$newname" != "$file" ]]; then
    echo "Renaming: $file → $newname"
    mv "$file" "$newname"
  fi
done

# in the ampseq environment

python /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/demux_nanopore_plates_edgesize20.py \
-p /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/kai_anopheles/plate_layout_kai.csv \
-b /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/kai_anopheles/internal_barcodes_kai.csv \
-f /mnt/storage11/sophie/kai/fastq \
-o samples.csv \
-m 0 \
-t 8

#3 Move the fastq files that have been double demultiplexed to another folder 
# and then continue with amplicon script.

for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done


#5 # Create a sample_file.csv
ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv 
#remove any unassigned/unclassified sample names here as you don't want to analyse those

# 6 some pre-run checks that are good to do
# Print any region in the bed file that does not match with an annotation in the gff file:
bedtools intersect -a AgamP4_chr.bed -b Anopheles_gambiae.AgamP4.56.gff3 -wa | sort | uniq > matched.bed
comm -23 <(cut -f1-4 AgamP4_chr.bed | sort) matched.bed
# this prints all regions that match between the bed file and the gff file to matched.bed and then prints any 
# that do not match in the terminal for you to see. You want the output to terminal to be blank, which indicates that all of the bed file regions are also in the gff.
# Check that regions in the bedfile match the reference fastq file
awk 'NR==FNR {chrlen[$1]=$2; next} $1 in chrlen && $3 <= chrlen[$1] {next} {print "Invalid:", $0}' Anopheles_gambiae.AgamP4.dna.toplevel.fa.fai AgamP4_chr.bed
# If nothing prints to terminal, you're good.


# 7 Make sure there is no samclip as this gets rid of all amplicon reads
# Mapping of amplicon data
python /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/sophie_nanopore_amplicon_script_minimap2_v2.py \
--index-file sample_file.csv \
--ref Anopheles_gambiae.AgamP4.dna.toplevel.fa \
--gff Anopheles_gambiae.AgamP4.56.gff3 \
--bed AgamP4_chr.bed \
--min-base-qual 30 \
--snpeff-db Anopheles_gambiae \
--threads 40 > amplicon_log.txt 2>&1


#8 Create amplicon x sample coverage matrix by combining coverage files (the sample_coverage_mean.txt files) to see how good they are across the sample set
