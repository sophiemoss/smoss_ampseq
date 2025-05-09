# The amplicon script from the group takes multiplexed files and processes them. 
# The following is designed to take files which you have already demultiplexed.
# This is useful if you are processing different samples at different times, and then want to run all together with the amplicon script for variant calling.

# Running amplicon script for demultiplexed fastq files. 

# You need to make an index file with the samples to process and no others, and use random barcode eg CCCCC because you have already demultiplexed these files

ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > anoph_gam_amp_index.txt

# Can manipulate in bash or export to pc and open in excel, add columns for sample, I1, I2, and CCCCC.
# Can save as CSV and import back to server to use.
# The amplicon script is designed for later annotation with snpeff, so it needs to use an ENSEMBL reference genome and gff file for this. gff file from https://metazoa.ensembl.org/.
# Amplicon script edited so that it would take in a samples list instead of using the index file.


# OPTION 1: call variants with both gatk and freebayes
# with mysamples.txt file, make this and take out Unknown and bu1033rs which are not needed/in the wrong format.
# amplicon_script5 is /mnt/storage8/sophie/amplicon-seq/scripts/amplicon_script5.py and in anopheles_ampseq repo

ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > mysamples.txt

python /mnt/storage8/sophie/amplicon-seq/scripts/amplicon_script5.py \
--read1 _1.fastq.gz \
--read2 _2.fastq.gz \
--ref /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
--gff /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3 \
--bed /mnt/storage8/sophie/bijagos_mosq_amplicon/AgamP4_chr.bed \
--position-info /mnt/storage8/sophie/bijagos_mosq_amplicon/agap.position_info.txt \
--samples-file mysamples.txt > amp_snpeff_log.txt 2>&1


# OPTION 2: using bam files produced, call variants with just freebayes
## Because this did not call all of the positions in the other samples, Jody ran the variant calling from the bam files that were produced using freebayes only:

freebayes -f /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
--haplotype-length -1 --pooled-continuous \
--targets /mnt/storage8/sophie/bijagos_mosq_amplicon/AgamP4_chr.bed *.bam > all_samples.vcf


# OPTION 3: call variants with gatk only
## Now doing with GATK only to see if there is a difference to freebayes:

# need to have all of the bam index files (.bam.bai) in the directory too
ls *bam | cat > bam.list

# GATK works by first using HaplotypeCaller to make genomic VCFs for each sample from their bam file, (genomic VCF does not contain genotype info)
# second using CombineGVCFs to combine these genomic VCFs into one file. Usually smaller file than the freebayes variant caller would produce.
# thirdly using GenotypeGVCFs to genotype the variants.


# Step 1) HaplotypeCaller
# -R is the reference.fasta, -L are the exome targets, -I are the input bam files, \
# -O is the output file -ERC GVCF more is used to add mutliple samples to the analysis
# can do this on multiple bam files at once by creating bam.list

gatk HaplotypeCaller \
-R /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
-L /mnt/storage8/sophie/bijagos_mosq_amplicon/AgamP4_chr.bed  \
-I bam.list \
-ERC GVCF \
-O bamfile.g.vcf \
> log.txt 2>&1

# now that we have the individual GVCFs make sure they are suffixed *.g.vcf if not specified in -O above

for i in *.vcf.idx; do mv -- "$i" "${i%.vcf.idx}.g.vcf.idx"; done

# but best practise for GATK is to run Haplotype Caller individually on all sampels using a loop, and then combine each of the GVCF files to a single GVCF file
# written a loop to go through each of the bam files and create a vcf for each of them using gatk haplotype caller
# flags for haplotype caller can be seen in the loop below.

ipython
import fastq2matrix as fm
import glob
from fastq2matrix import run_cmd

for i in glob.glob('./*.bam'):

    Joy = i.replace('.bam','.vcf')
    
    run_cmd (f"gatk HaplotypeCaller -R /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
    -L /mnt/storage8/sophie/bijagos_mosq_amplicon/AgamP4_chr.bed -I {i} -ERC GVCF -O {Joy}")

# I want it to include all sites so that I can see in the filter column why a site was not included previously
#--output-mode EMIT_ALL_CONDIFENT_SITES which produces calls at variant sites and confident reference sites

for i in glob.glob('./*.bam'):
    
    Joy = i.replace('.bam','.emit.vcf')
    
    run_cmd (f"gatk HaplotypeCaller -R /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
    -L /mnt/storage8/sophie/bijagos_mosq_amplicon/AgamP4_chr.bed -I {i} -ERC GVCF --output-mode EMIT_ALL_CONFIDENT_SITES -O {Joy} > log.txt 2>&1")


# Step 2) now we use CombineGVCFs to produce a combined multi-sample gVCF. 

 gatk CombineGVCFs \
   -R /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
   --variant input.list \
   -O combined.g.vcf.gz

# Step 3) then use GenotypeGVCFs to perfom joint genotyping on all of the samples

 gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /mnt/storage8/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
   -V combined.g.vcf.gz \
   --include-non-variant-sites \
   -O combined_genotype_gatk_INVS.vcf.gz > log.txt 2>&1

   # Known issues:
   # When calling variants, freebayes correctly keeps in the denominator, maintaing the homozygous reference calls for samples where a SNP is not present
   # GATK alone and GATK+Freebayes combined do not do this, the homozygous reference calls are filtered out at the HaplotypeCaller stage
  



