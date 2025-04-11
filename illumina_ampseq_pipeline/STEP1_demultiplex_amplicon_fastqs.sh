# Demultiplexing multiple fastq files using the fastq2matrix demultiplex.py script

# You need the index files which include data on which sample is tagged with which barcode
# Moving multiple index files at once
cat numbers.txt | parallel "cp /mnt/storage8/sophie/bijagos_mosq_amplicon/ampseq_output/anpoolseq{}_output/anpoolseq{}_index.csv ." > indexlog.txt 2>%1

# Run the demultplex script

demultiplex_fastq.py --R1 %(read1)s --R2 %(read2)s --index %(index_file)s

# For a single sample

demultiplex_fastq.py --R1 Anpoolseq1_R1_001.fastq.gz --R2 Anpoolseq1_R2_001.fastq.gz --index anpoolseq1_index.csv > log.txt 2>&1

# For parallel samples

ls *R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' > samples.txt

cat samples.txt | parallel -j 2 "demultiplex_fastq.py --R1 {}_R1_001.fastq.gz --R2 {}_R2_001.fastq.gz --index {}_index.csv" > log.txt 2>&1

# Once all of the raw fastq files are demultiplexed and in the same directory, you can run the amplicon pipeline on them all at once.

