#5 

ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv

# 6 some pre-run checks that are good to do

# Print any region in the bed file that does not match with an annotation in the gff file:
bedtools intersect -a AgamP4_chr.bed -b Anopheles_gambiae.AgamP4.56.gff3 -wa | sort | uniq > matched.bed
comm -23 <(cut -f1-4 AgamP4_chr.bed | sort) matched.bed
# this prints all regions that match between the bed file and the gff file to matched.bed and then prints any 
# that do not match in the terminal for you to see. You want the output to terminal to be blank, which indicates that all of the bed file regions are also in the gff.

# Check that regions in the bedfile match the reference fastq file
awk 'NR==FNR {chrlen[$1]=$2; next} $1 in chrlen && $3 <= chrlen[$1] {next} {print "Invalid:", $0}' Anopheles_gambiae.AgamP4.dna.toplevel.fa.fai AgamP4_chr.bed
# If nothing prints to terminal, you're good.

#6 
# Make sure there is no samclip as this gets rid of all amplicon reads
# Mapping of amplicon data

python /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/sophie_nanopore_amplicon_script_minimap2_v2.py \
--index-file sample_file.csv \
--ref AnoDarl_H01.genomic.fasta \
--gff GCF_943734745.1_idAnoDarlMG_H_01_genomic.gff \
--bed darlingi.bed \
--min-base-qual 30 \
--threads 40 \
--snpeff-db Anopheles_darlingi > amplicon_log.txt 2>&1 

#Anopheles_darlingi_2 snpeff database I ran separately for Aminul - but this is not available in the ampseq environment
# Now just checking that it all runs in the ampseq environment for future work (with the script that includes snpeff)

#7 Create amplicon x sample coverage matrix by combining coverage files (the sample_coverage_mean.txt files) to see how good they are across the sample set



