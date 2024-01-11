# FILTERING OF COMBINED GENOTYPED VCF
# same filters used for both gambiae pipeline and melas pipeline

########### FILTER 1: remove INDELS with bcftools, here -M2 indicates that bcftools should split multi-allelic sites into multiple biallelic sites, 
# keeping this information
# -m2 is used in conjunction with -M2 to apply the minor-allele-based decomposition. 
# This means that bcftools will decompose multi-allelic sites using the minor allele as the reference allele in the biallelic split.
# here -v tells bcftools to only view SNPS, so indels are excluded

# bcftools view -M2 -m2 -v snps melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz -Oz -o bi_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz
tabix -p vcf bi_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

echo "indexed vcf after bcftools view -M2 -m2 -v snps to make bi_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz.tbi"

# filter out the contigs from the VCF file, note that to produce a properly bgzipped vcf file you need the -Oz flag
bcftools view bi_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz --regions 2L,2R,3L,3R,Mt,X,Y_unplaced | bcftools sort -Oz -o bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

echo "filtered out contigs from VCF"

# view the unique chromosome names
bcftools query -f '%CHROM\n' bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt
tabix -p vcf bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

echo " viewed unique chromosome names and printed them to unique_chromosomes_filtered.txt"
echo "indexed vcf which has had the contigs filtered out"

########### FILTER 2: Filter samples to keep those with 40% of genome with > 10x coverage, and min-ac=1 so that all variants that remain are still variants after sample removal
# I have identified samples that are above this threshold using my basic WGS stats
# create file with samples to keep: sample_40_10.txt
bcftools view -S samples_40_10.txt bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz --min-ac=1 -Oz -o miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz
tabix -p vcf miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

## count variants, gambiae: 21679588
line_count=$(bcftools view -H miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | wc -l)

echo "The number of variants after including only samples with 40% 10x coverage across genome, and --min-ac=1 is: $line_count"

########### FILTER 3: for GATK standard filtering
# GATK filter recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# excellent explanation of this hard filtering and additional soft filtering
# https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/#ex3.4

echo "Tagging variants for gatk Variant Filtration"

gatk VariantFiltration \
-R /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
-V miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz \
-filter "QD < 5.0" --filter-name "QD5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_tagged_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

## note there is a repeated warning here JEXL engine undefined variable, not a problem
## because this is coming from where there are positions with no coverage
## https://gatk.broadinstitute.org/hc/en-us/community/posts/4408733963803-GATK-Variant-Filtration-undefined-variable

## remove the SNPs that have been tagged with the above filters
bcftools view -f 'PASS' gatk_tagged_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz -Oz -o gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz
tabix -p vcf gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

echo "Removed variants that did not meet gatk filtering thresholds"

## count variants, gambiae: 20037931
line_count2=$(bcftools view -H gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | wc -l)

echo "The number of variants is: $line_count2"
echo "Hard filering complete, now moving to soft filtering"

## ADDITIONAL SOFT FILTERING

## After hard quality filtering, we have a data set containing only variant sites that we trust with high confidence. 
## However, so far we have only removed variants that had skewed values across all samples – 
## we haven’t yet looked at the individual genotypes at each variant site. 
## For example, we may have kept a variant that has a high quality across all of the individual samples, but there 
## may be some individuals where the read depth (FMT/DP) or the genotype quality (GQ) for the individual genotype is very low, 
## so we don’t trust the calls for these individuals

########### FILTER 4: Filter reads that have a read depth below 5 OR a genotype quality below 20
## The difference is, that instead of the whole variant site, we are now considering single genotypes 
## (the information in the “FORMAT” fields”) using -e 'FMT/DP<3 | FMT/GQ<20' and we will not remove the whole variant site, 
## but only set the respective genotypes to missing (./.) by using the bcftools filter -S . command.
## Here we use the logical operator “|” instead of “||”, since using “||” would mean that every genotype at a variant site 
## is set to missing, even if only one genotype doesn’t fulfill the threshold
bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20' gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz -Oz -o DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

## check this, the genotype here should have been set to missing ./. as we are filtering for depth below 5
#bcftools query -i 'FMT/DP<5' -f '[GT=%GT\tDP=%DP\tGQ=%GQ\t]\n' DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | less -S

## happy? index vcf
tabix -p vcf DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

echo "Filtered for FMT/DP<5 and FMT/GQ<20"

## count variants, gambiae: 20037931. Same number as all are filters above just set to missing, so genotypes removed but site still in VCF.
line_count3=$(bcftools view -H DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | wc -l)
echo "The number of variants is $line_count3"

########### FILTER 5: exclude -e all sites at which no alternative alleles are called for any of the samples
bcftools filter -e 'AC==0' DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz -O z -o AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

tabix -p vcf AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

echo "Excluded all sites where there are no alternative alleles called anymore due to missingness"

## count variants
line_count4=$(bcftools view -H AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | wc -l)
echo "The number of variants after AC==0 filter is $line_count4"

########### FILTER 6
## Remove variants with a high amount of missing genotypes and filter on minor allele frequency
## The data set overall will now have lots of missing data, because we have replaced calls with depth <5 or quality below 20 with ./.  
## Therefore we will now remove all variants that have more than 20% missing genotypes or MAF < 0.01
## Filtering for MAF 0.01 means that remaining sites in the VCF need to have a minimum minor allele frequency of 1%
## This dataset includes 44 individuals (thus maximally 88 alleles per site). If we filter for MAF 0.01, 
## a variant is excluded if its minor allele is found in less than 0.88 genotypes (1% of 88). 
## This would remove alleles where the minor allele occurs 0.88 times or less, which for this particular dataset means a 
## variant would be removed if there are no minor alleles? So if there is no alternate allele the variant is removed.

bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz -O z -o FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

line_count5=$(bcftools view -H FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | wc -l)
echo "The number of variants after AC==0 filter is $line_count5"
echo "Filtering complete, variants went from $line_count to $line_count2 to $line_count3 to $line_count4 to $line_count5"

# Filtering complete!

########## PHASE VCF FILE ###########

# Phase the filtered vcf using beagle, https://faculty.washington.edu/browning/beagle/beagle_5.2_13Oct21.pdf
# necessary for some selection statistical analyses
# could also use phasing pipeline from malariagen. Have a look at this, and understand gatk parameters above.

# use beagle conda environment. mamba install beagle.
conda activate beagle
beagle -Xmx500g gt=FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz out=2019melasvcfphased

## Check the number of SNPs in the phased and unphased VCF files, these should be the same

#OLD
# unphased melas
bcftools query -f '%CHROM\t%POS\n' FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz | awk '$1=="2L"' | wc -l
# 2652491

## phased
# bcftools query -f '%CHROM\t%POS\n' 2019melasvcfphased.vcf.gz | awk '$1=="2L"' | wc -l
# 2652491
