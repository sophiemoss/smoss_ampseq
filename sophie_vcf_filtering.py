#script to take combined vcf file and look for SNPs, INDELs across sample.
#includes filtering section
#this is saved in /mnt/storage8/sophie/amplicon-seq/scripts


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

## Filtering steps

# Import necessary python packages
import pandas as pd
import numpy as np

### combined SNPs
# Read in the dataframe as VarDf and then add column headers

VarDf = pd.read_csv('combined_snps_trans.txt',header=None,sep='\t')
VarDf.columns = 'SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\tANN'.split('\t')

#check if all calls are in VarDf
#eg. making sure the variant caller has kept in positions across all samples where one sample has a SNP
VarDf.to_csv('python_all_snps_check1.tsv',index=None,sep='\t')

# Filter out the rows where is no call
VarDfFilt = VarDf[VarDf['GT'].isin(['./.','.'])==False]

#check if all calls are in VarDfFilt
VarDfFilt.to_csv('python_snps_check2.tsv',index=None,sep='\t')

# Extract Allele Depths = Remember GATK removes uninformative reads which causes difference in SUM(AD) vs DP and also '.'==0
# change these commands to use .loc

VarDfFilt['REF_Depth'] = [int(X.replace('.','0').split(',')[0]) for X in VarDfFilt['AD'].tolist()]

VarDfFilt['ALT_Depth'] = [int(X.replace('.','0').split(',')[1]) for X in VarDfFilt['AD'].tolist()]

# Calculate Sum of (AD) columns. Note replacing '.' with 0. 
VarDfFilt['Sum of AD'] = VarDfFilt['REF_Depth'] + VarDfFilt['ALT_Depth']

# Converts DP values to integers (if you want to compare)
VarDfFilt['DP'] = [int(X) for X in VarDfFilt['DP']]

# To filter to keep only the variants with matching DP and Sum of AD. (DECIDE IF UNINFORMATIVE READS ARE IMPORTANT)
# VarDfFilt = VarDfFilt[VarDfFilt['Sum of AD']==VarDfFilt['DP']]

# Standardise GT Column
VarDfFilt['GT'] = VarDfFilt['GT'].replace({'1/1':1,
                                            '0/1':0.5,
                                            '1/0':0.5,
                                            '0/0':0})
#Create CSV with SNPs so far
VarDfFilt.to_csv('python_combined_snps_alldepths.tsv',index=None,sep='\t')

# Filter out the rows where DP is less than 30, leaving only DP >= 30
VarDfFilt2 = VarDfFilt[VarDfFilt['DP']>=30]

# Make a separate dataframe of the SNPs with DP less than 30 just in case needed
VarDPUnder30 = VarDfFilt[VarDfFilt['DP']<30]

# Write the dateframe out as tsv - do not do a CSV because there are comments which will get confusing
VarDfFilt2.to_csv('python_combined_snps_removednocallsand<30_.tsv',index=None,sep='\t')

# Test if multiple values in the ALT column seperated by a comma
X = VarDfFilt[VarDfFilt['ALT'].str.contains(',')]

## add genotype based on 25% and 75% thresholds for heterozygous
    ## add reference depth as % of total depth column, and alternate depth as % of total depth as column
VarDfFilt2['RefProp'] = VarDfFilt2["REF_Depth"]/VarDfFilt2["Sum of AD"]
VarDfFilt2['AltProp'] = VarDfFilt2["ALT_Depth"]/VarDfFilt2["Sum of AD"]
    ## if REF depth is 75% or more of 'Sum of AD', call is homozygous reference
    ## if ALT depth is 75% or more of 'Sum of AD', call is homozygous alternate
    ## if ALT depth is between 25% and 75%, then call is heterozygous

# add a new empty column at column position 17
VarDfFilt2.insert(17,"GT_new","")

# write a for loop to add new genotype call
# print Mistake? if there is a probable error in the GT_new column
for index, row in VarDfFilt2.iterrows():
    if row['RefProp'] >= 0.75:
        VarDfFilt2.at[index, 'GT_new'] = 'Homo_Ref'
    elif row['AltProp'] >= 0.75:
        VarDfFilt2.at[index, 'GT_new'] = 'Homo_Alt'
    elif 0.25 <= row['AltProp'] < 0.75:
        VarDfFilt2.at[index, 'GT_new'] = 'Hetero'
    elif row['RefProp'] or ['AltProp'] == 'NaN':
        VarDfFilt2.at[index, 'GT_new'] = 'NoCall'
    else:
        VarDfFilt2.at[index, 'GT_new'] = 'Mistake?'

#check the GT_new calls
VarDfFilt2['GT_new'] = 'NA'

#If there is a Mistake? then print these rows to look more closely
#other_rows = VarDfFilt2[VarDfFilt2['GT_new'] == 'Mistake?']
#print(other_rows)

###ALERT###
#the NaN rows are where there is no call - this is VERY important: MUST use the GT_new column. The old GT column has 0.0 for positions which are homo ref AND empty. GT_new solves this.

## add gene name as separate column (3rd column in ANN, remember counts from 0)
VarDfFilt2['Gene_name'] = [(X.split('|')[3]) for X in VarDfFilt2['ANN'].tolist()]

## add gene ID as a separate column (4th column in ANN, remember counts from 0)
VarDfFilt2['GeneID'] = [(X.split('|')[4]) for X in VarDfFilt2['ANN'].tolist()]

## write CSV with GT_new
VarDfFilt2.to_csv('python_combined_snps_filtered_genotyped.tsv',index=None,sep='\t')

## annotate specific variants of interest in a separate column

# read in specific genomic locations of interest
snpsofinterest = pd.read_csv('anopheles_gambiae_sl_snps.tsv',sep='\t')
# match position from snpsofinterest with POS column in VarDfFilt2
# add empty column to VarDfFilt2 to hold the values
VarDfFilt2['snpofinterest'] = 'NA'
# iterate through each row in VarDfFilt2
for index, row in VarDfFilt2.iterrows():
    # check if the 'POS' value in VarDfFilt2 matches a 'position' value in snpsofinterest
    snp_match = snpsofinterest[snpsofinterest['position'] == row['POS']]
        # if a match is found, set the 'snpofinterest' value in VarDfFilt2
    if not snp_match.empty:
        var_value = f"{snp_match['position'].iloc[0]}_{snp_match['gene'].iloc[0]}_{snp_match['AA.change'].iloc[0]}"
        VarDfFilt2.at[index, 'snpofinterest'] = var_value

## add species and island of samples using metadata

## read in metadata with the geographic information
metadata_anopheles_2019 = pd.read_csv('2019_anopheles_metadata.tsv',sep='\t')

## add empty columns to VarDfFilt2 to hold the species and island values
VarDfFilt2['species'] = 'NA'
VarDfFilt2['island'] = 'NA'

# change SAMPLE names in VarDFFilt2 to all be lowercase letters
VarDfFilt2['SAMPLE'] = VarDfFilt2['SAMPLE'].str.lower()

# iterate through each row in VarDfFilt2
for index, row in VarDfFilt2.iterrows():
    # check if the 'SAMPLE' value in VarDfFilt2 matches a 'Sample' value in metadata_anopheles_2019
    sample_match = metadata_anopheles_2019[metadata_anopheles_2019['Sample'] == row['SAMPLE']]
        # if a match is found, set the 'species' value in VarDfFilt2
    if not sample_match.empty:
        var_value = f"{sample_match['Species'].iloc[0]}"
        VarDfFilt2.at[index, 'species'] = var_value

 # iterate through each row in VarDfFilt2
for index, row in VarDfFilt2.iterrows():
    # check if the 'SAMPLE' value in VarDfFilt2 matches a 'Sample' value in metadata_anopheles_2019
    sample_match = metadata_anopheles_2019[metadata_anopheles_2019['Sample'] == row['SAMPLE']]
        # if a match is found, set the 'island' value in VarDfFilt2
    if not sample_match.empty:
        var_value = f"{sample_match['Island'].iloc[0]}"
        VarDfFilt2.at[index, 'island'] = var_value

# save dataframe as tsv
VarDfFilt2.to_csv('python_combined_snps_filtered_genotyped_labelled.tsv',index=None,sep='\t')
VarDfFilt2 = pd.read_csv('python_combined_snps_filtered_genotyped_labelled.tsv',sep='\t')

#filter for quality of snps >30
VarDfFilt2 = VarDfFilt2[VarDfFilt2['QUAL']>=30]
VarDfFilt2.to_csv('python_combined_snps_filtered_genotyped_labelled_qualfiltered.tsv',index=None,sep='\t')


### Combined INDELS

# Read in the dataframe as InDf and then add column headers

InDf = pd.read_csv('combined_indels_trans.txt',header=None,sep='\t')
InDf.columns = 'SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\tANN'.split('\t')