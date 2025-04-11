# This script takes a combined vcf file from the output of the amplicons sequencing pipeline
# and generates a SNPs and an INDELS VCF, which is annotated using snpeff
# the script then pulls the SNPS, INDELS and annotation out into txt files

#import necessary python packages

import fastq2matrix as fm
from fastq2matrix import run_cmd

run_cmd(r"bcftools query -f '%CHROM\t%POS[\t%DP]\n' combined_genotype_gatk.vcf.gz > tmp.txt")

#FMT/DP>10 means that both FMT and DP conditions need to be satisfied. Depth must be > 10 and FMT means format
run_cmd("bcftools filter -i 'FMT/DP>10' -S . combined_genotype_gatk.vcf.gz | bcftools sort -T . | bcftools norm -m - -Oz -o tmp.vcf.gz")
run_cmd("bcftools view -v snps tmp.vcf.gz > snps.vcf")
run_cmd("bgzip snps.vcf")
run_cmd("snpEff Anopheles_gambiae snps.vcf.gz > snps.ann.vcf")
run_cmd("bgzip snps.ann.vcf")
run_cmd("tabix snps.vcf.gz")
run_cmd("tabix snps.ann.vcf.gz")
run_cmd("bcftools view -v indels tmp.vcf.gz > indels.vcf")
run_cmd("bgzip indels.vcf")
run_cmd("snpEff Anopheles_gambiae indels.vcf.gz > indels.ann.vcf")
run_cmd("bgzip indels.ann.vcf")
run_cmd("tabix indels.vcf.gz")
run_cmd("tabix indels.ann.vcf.gz")

run_cmd(r"bcftools query snps.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%AD\t%DP\n]' > combined_snps.txt")    
run_cmd(r"bcftools query snps.ann.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%AD\t%DP\t%ANN\n]' > combined_snps_trans.txt")
run_cmd(r"bcftools query indels.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%AD\t%DP\n]' > combined_indels.txt")
run_cmd(r"bcftools query indels.ann.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%AD\t%DP\t%ANN\n]' > combined_indels_trans.txt")
