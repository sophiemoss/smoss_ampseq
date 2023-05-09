# The amplicon pipeline results in a combined VCF file of all samples.
# Txt files of SNPs and indels are then made from this combined VCF.
# This script has been developed to filter the text files of SNPs and INDELs for downstream processing
# parameters and metadata will need to be changed by the user depending on their dataset

# You need to install pandas. You can do this with conda install -c anaconda pandas

## Filtering steps

# Import necessary python packages
import pandas as pd
import numpy as np

### Combined SNPs
# Read in the dataframe as VarDf and then add column headers

VarDf = pd.read_csv('combined_snps_trans.txt',header=None,sep='\t')
VarDf.columns = 'SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\tANN'.split('\t')

#check if all calls are in VarDf
#eg. making sure the variant caller has kept in positions across all samples where one sample has a SNP
VarDf.to_csv('1_python_all_snps_check.tsv',index=None,sep='\t')

# Filter out the rows where is no call
VarDfFilt = VarDf[VarDf['GT'].isin(['./.','.'])==False]

#check if all calls are in VarDfFilt
VarDfFilt.to_csv('2_python_snps_check2.tsv',index=None,sep='\t')

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
VarDfFilt.to_csv('3_python_combined_snps_alldepths.tsv',index=None,sep='\t')

# Filter out the rows where DP is less than 30, leaving only DP >= 30
VarDfFilt2 = VarDfFilt[VarDfFilt['DP']>=30]

# Make a separate dataframe of the SNPs with DP less than 30 just in case needed
VarDPUnder30 = VarDfFilt[VarDfFilt['DP']<30]

# Write the dateframe out as tsv - do not do a CSV because there are comments which will get confusing
VarDfFilt2.to_csv('4_python_combined_snps_removednocallsandDP<30_.tsv',index=None,sep='\t')

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
VarDfFilt2['GT_new'] = 'NA'

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
VarDfFilt2.to_csv('5_python_combined_snps_filtered_0.250.75genotyped.tsv',index=None,sep='\t')

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
VarDfFilt2.to_csv('6_python_combined_snps_filtered_genotyped_labelled.tsv',index=None,sep='\t')
VarDfFilt2 = pd.read_csv('6_python_combined_snps_filtered_genotyped_labelled.tsv',sep='\t')

#filter for quality of snps >30
VarDfFilt2 = VarDfFilt2[VarDfFilt2['QUAL']>=30]
VarDfFilt2.to_csv('7_python_combined_snps_filtered_genotyped_labelled_qualfiltered.tsv',index=None,sep='\t')


### Now move on to filtering the combined INDELS.

# Read in the dataframe as InDf and then add column headers

# InDf = pd.read_csv('combined_indels_trans.txt',header=None,sep='\t')
# InDf.columns = 'SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\tANN'.split('\t')