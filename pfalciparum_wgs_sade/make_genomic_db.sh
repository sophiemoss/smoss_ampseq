## to make a genomics database of sample VCFs, use the following
ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > fastq2vcfsamples.txt

# Create DB and import new isolates to store together information from multiple vcfs
# also creates .json file

python /mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import \
--sample-file fastq2vcfsamples.txt \
--ref Pfalciparum.genome.fasta \
--prefix sade_falciparum_2023 \
--vcf-dir .\
> mergevcf_import_log.txt 2>&1

## now merge VCF files

python /mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py genotype\
 --ref Pfalciparum.genome.fasta \
 --prefix sade_falciparum_2023 \
  > mergevcf_genotype_log.txt 2>&1

# resulting vcf is called sade_falciparum_2023 genotyped vcf.gz