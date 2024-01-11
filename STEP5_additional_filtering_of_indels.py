### Now move on to filtering the combined INDELS.

import pandas as pd
import numpy as np

# Read in the dataframe as InDf and then add column headers

InDf = pd.read_csv('combined_indels_trans.txt',header=None,sep='\t')
InDf.columns = 'SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\tANN'.split('\t')

# Filter out where there is no depth
InDf = InDf[InDf['DP'].isin(['./.','.'])==False]

# Converts DP values to integers (if you want to compare)
InDf['DP'] = [int(X) for X in InDf['DP']]

# Filter out the rows where DP is less than 30, leaving only DP >= 30
InDf2 = InDf[InDf['DP']>=30]

#filter for quality of calls >30
InDf2 = InDf2[InDf2['QUAL']>=30]

# Write to csv

InDf2.to_csv('8_indels_depth_filtered.tsv',index=None,sep='\t')