
######################## SELECTION STATISTICS #########################

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps')
os.getcwd()

## A.miles
# IHS scans were
# computed following methods described in Voight et al. (2006) as implemented in scikit-allel
# version 0.21.1. For each population, SNPs with minor allele frequency above 5% were used.
# IHS scores were normalised within each chromosome (2, 3, X).

# %%
## have now phased VCFs using beagle

import numpy as np
import allel
import zarr
import pandas as pd
import sys

## convert phased, filtered, VCF file to zarr file
# %%
# allel.vcf_to_zarr('2022gambiaevcfphased.vcf.gz', '2022gambiaevcfphased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('2022gambiaevcfphased.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

# %%
## import metadata
df_samples=pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count

# %% select resistant samples from metadata by index value and store as array

res_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
res_samples

# %% select genotypes for variants for those resistant samples and store as a genotype dask array

gt_res_samples = gt.take(res_samples, axis=1)
gt_res_samples

# %% select variants that are segregating within res_samples as only these will be informative
## also some selection tests don't support multiallelic variants, so just keep biallelics
## for this pipeline the VCF is already filtered so should be no biallelic SNPs anyway

# %% compute allele counts for resistant samples and create allele counts array
ac_res = gt_res_samples.count_alleles(max_allele=8).compute()
# %% filter for those that are segregating and biallelic and store as a boolean array
# is.segregating() just finds variants where more than one allele is observed. is_non_segregating() finds non-segregating variants (where at most one allele is observed)

res_seg_variants = ac_res.is_segregating() & ac_res.is_biallelic_01()
# %% remove variants that are on Y_unplaced using a boolean mask
chrom = callset['variants/CHROM'][:]
exclude_chrom = 'Y_unplaced'
res_seg_variants = res_seg_variants & (chrom != exclude_chrom)
# %% make an allele counts array from this
ac_res_seg = ac_res.compress(res_seg_variants, axis=0)
# %% also make a genotype dask array
gt_res_seg = gt_res_samples.compress(res_seg_variants, axis = 0)
gt_res_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_res_seg = gt_res_seg.to_haplotypes().compute()
h_res_seg

# %% we need variant positions
pos = callset['variants/POS'][:]
pos_res_seg = pos.compress(res_seg_variants, axis=0)
pos_res_seg

# %% also store chromosome of each variant as we need this for shading the plots later
chrom = callset['variants/CHROM'][:]
chrom_res_seg = chrom.compress(res_seg_variants, axis=0)
chrom_res_seg

# %%
# some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_res_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %%
# compute raw iHS

ihs_res_raw = allel.ihs(h_res_seg, pos_res_seg, use_threads=True, include_edges=True)
ihs_res_raw

# %%
%matplotlib inline
import matplotlib.pyplot as plt
from datetime import datetime

# %% Plot raw iHS

fig, ax = plt.subplots()
ax.hist(ihs_res_raw[~np.isnan(ihs_res_raw)], bins=20)
ax.set_xlabel('Raw IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize iHS. Note that iHS has been calculated with unpolarized data, so only the magnitude of iHS

ihs_res_std = allel.standardize_by_allele_count(ihs_res_raw, ac_res_seg[:, 1])
ihs_res_std

# %% Generate timestamp with current date and time for the figure you are about to make
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

# %% Plot standardized iHS over genome and display which iHS values are significant
# Default red line for significance for iHS is 4. Default for XP-EHH is 5. Defauly for RSB is 5.
# I am using iHS significanc of 5 because otherwise it is too many SNPs at 4.

#fig, ax = plt.subplots(figsize=(10, 3))
#ax.plot(pos_res_seg, np.abs(ihs_res_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
#ax.axhline(y=5, color='red', linestyle='--')
#ax.set_xlabel('Genomic position (bp)')
#ax.set_ylabel('$|IHS|$')
#ax.set_ylim(0, 9)
#ax.legend()
#

# %% Plot standardized iHS across genome and shade chromosomes

# Shading manhattan based on chromosome
chromosome_colours = {
    '2L': 'red',
    '2R': 'orange',
    '3L': 'yellow',
    '3R': 'green',
    'X': 'blue',
    'Mt': 'black',
    'Y_unplaced': 'purple',
}

fig, ax = plt.subplots(figsize=(10, 3))

offset = {
    '2L':0,
    '2R':49364325,
    '3L':110909430,
    '3R':152872865,
    'X':206073549,
    'Mt':230466657,
    'Y_unplaced':230482020
}

# Create a list of chromosome names and their corresponding positions
chromosome_labels = [
    ('2L', 25000000),
    ('2R', 75000000),
    ('3L', 130000000),  # Adjust as needed
    ('3R', 180000000),  # Adjust as needed
    ('X', 218000000),   # Adjust as needed
    ('Mt', 230466657),
    ('Y_unplaced', 230482020)
]

# Extract chromosome names and positions from the list
chromosomes, positions = zip(*chromosome_labels)

# Set the x-axis ticks and labels
ax.set_xticks(positions)
ax.set_xticklabels(chromosomes)

# Loop through chromosomes and plot variants with different colors
for chrom, color in chromosome_colours.items():
    is_chrom = (chrom_res_seg == chrom)
    xpos = offset[chrom] + pos_res_seg[is_chrom]
    ax.plot(xpos, np.abs(ihs_res_std[0][is_chrom]), linestyle=' ', marker='o', mfc='none', mew=.5, mec=color, label=chrom)

# Plot variants above the significant line in different colour
above_significance = np.abs(ihs_res_std[0]) > 5  # You can adjust the threshold here
#ax.plot(xpos[above_significance], np.abs(ihs_res_std[0][above_significance]), linestyle=' ', marker='o', mfc='cyan', mew=0.25, mec='blue', label='Above Significance')

ax.axhline(y=5, color='red', linestyle='--')
ax.set_xlabel('Chromosome')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9)
ax.legend()

# Disable scientific notation for x-axis so that full numbers are printed
#ax.get_xaxis().get_major_formatter().set_scientific(False)

# Save figure
#timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
#filename = f'shaded_standardised_ihs_histogram_with_cutoff_{timestamp}.png'
#plt.savefig(filename)

# %% Find the index of the variant with the highest iHS value
idx_hit_max = np.nanargmax(ihs_res_std[0])

# %% Find the genomic position of top hit
pos_res_seg[idx_hit_max]
print(f'Genomic position with highest iHS value:', pos_res_seg[idx_hit_max])

# %% Visualise EHH decay around top hit with a Voight plot
# pull out haplotypes for region around top hit
flank_size = 2000
h_hit = h_res_seg[idx_hit_max - flank_size:idx_hit_max + flank_size + 1]
h_hit

# %%
fig = allel.fig_voight_painting(h_hit[:, h_res_seg[idx_hit_max] == 0], index=flank_size, height_factor=0.02)
fig.suptitle('Reference allele', y=1);

# %% 
fig = allel.fig_voight_painting(h_hit[:, h_res_seg[idx_hit_max] == 1], index=flank_size, height_factor=0.02)
fig.suptitle('Alternate allele', y=1);

# %% list all positions with iHS value over certain threshold and save to a file

res_ihs_positions_above_threshold_6 = pos_res_seg[ihs_res_std[0] >= 6]
res_ihs_positions_above_threshold_6

with open("res_ihs_positions_above_threshold_6.txt", "w") as file:
    for position in res_ihs_positions_above_threshold_6:
        file.write(str(position) + "\n")
# %%
res_ihs_positions_above_threshold_5 = pos_res_seg[ihs_res_std[0] >= 5]
res_ihs_positions_above_threshold_5

with open("res_ihs_positions_above_threshold_5.txt", "w") as file:
    for position in res_ihs_positions_above_threshold_5:
        file.write(str(position) + "\n")

# %% ##################### SUSCEPTIBLE POPULATION #################

### Now analyse iHS in susceptible population

sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
sus_samples

# %%
## select genotypes for resistant samples

gt_sus_samples = gt.take(sus_samples, axis=1)
gt_sus_samples

# %%
## select variants that are segregating within sus_samples as only these will be informative

ac_sus = gt_sus_samples.count_alleles(max_allele=8).compute()
sus_seg_variants = ac_sus.is_segregating() & ac_sus.is_biallelic_01()
ac_sus_seg = ac_sus.compress(sus_seg_variants, axis=0)
gt_sus_seg = gt_sus_samples.compress(sus_seg_variants, axis = 0)
gt_sus_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_sus_seg = gt_sus_seg.to_haplotypes().compute()
h_sus_seg

# %%
## compute iHS values

# we need variant positions
pos = callset['variants/POS'][:]
pos_sus_seg = pos.compress(sus_seg_variants, axis=0)
pos_sus_seg

# %%
# some variants in 1000 genomes have multiple variants at the same genomic position, which causes problems for some selection tests in scikit-allel. Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_sus_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %%
# compute raw iHS

ihs_sus_raw = allel.ihs(h_sus_seg, pos_sus_seg, use_threads=True, include_edges=True)
ihs_sus_raw

# %%

%matplotlib inline
import matplotlib.pyplot as plt

# %%

fig, ax = plt.subplots()
ax.hist(ihs_sus_raw[~np.isnan(ihs_sus_raw)], bins=20)
ax.set_xlabel('Raw IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize iHS

ihs_sus_std = allel.standardize_by_allele_count(ihs_sus_raw, ac_sus_seg[:, 1])

# %% 
fig, ax = plt.subplots()
ax.hist(ihs_sus_std[0][~np.isnan(ihs_sus_std[0])], bins=20)
ax.set_xlabel('Standardized IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% note that iHS has been calculated with unpolarized data, so only the magnitude of iHS
# is informative, not the sign.
ihs_sus_std

# %% plot over the genome
fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_sus_seg, np.abs(ihs_sus_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9);

# %% look for the index of the variant with the highest iHS signal
idx_hit_max = np.nanargmax(ihs_sus_std[0])
idx_hit_max

# %% genomic position of top hit
pos_sus_seg[idx_hit_max]

# to check for variant in VCF use: bcftools view 2022gambiaevcfphased.vcf.gz | grep '28467703'

# %% Visualise EHH decay around top hit with a Voight plot
# pull out haplotypes for region around top hit
flank_size = 2000
h_hit = h_sus_seg[idx_hit_max - flank_size:idx_hit_max + flank_size + 1]
h_hit

# %%
fig = allel.fig_voight_painting(h_hit[:, h_sus_seg[idx_hit_max] == 0], index=flank_size, height_factor=0.02)
fig.suptitle('Reference allele', y=1);

# %% 
fig = allel.fig_voight_painting(h_hit[:, h_sus_seg[idx_hit_max] == 1], index=flank_size, height_factor=0.02)
fig.suptitle('Alternate allele', y=1);

# %%
# Plot iHS values and include line of which iHS values are significant

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_sus_seg, np.abs(ihs_sus_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.axhline(y=5, color='red', linestyle='--')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9)
ax.legend()

# %% list all positions with iHS value over certain threshold

sus_ihs_positions_above_threshold_6 = pos_sus_seg[ihs_sus_std[0] >= 6]
with open("sus_ihs_positions_above_threshold_6.txt", "w") as file:
    for position in sus_ihs_positions_above_threshold_6:
        file.write(str(position) + "\n")

sus_ihs_positions_above_threshold_5 = pos_sus_seg[ihs_sus_std[0] >= 5]
with open("sus_ihs_positions_above_threshold_5.txt", "w") as file:
    for position in sus_ihs_positions_above_threshold_5:
        file.write(str(position) + "\n")

# %% ########### Cross-population extended haplotype homozygosity (XPEHH) ###########

## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
# create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS
## VCF is phased so we can convert genotype arrays made earlier to haplotype array

## A.miles
# XPEHH scans were computed following methods described in Sabeti et al. (2007)
# as implemented in scikit-allel version 0.21.1. For each population comparison, SNPs
# with a minor allele frequency greater than 5% in the union of both populations were
# used. XPEHH scores were normalised within each chromosome (2, 3, X).

h_sus = gt_sus_samples.to_haplotypes().compute()
h_sus

h_res = gt_res_samples.to_haplotypes().compute()
h_res

ac_gt = gt.count_alleles(max_allele=8).compute()

# %%
# get variant positions

pos = callset['variants/POS'][:]
pos

# %% look at shapes of arrays
print("h_sus shape:", h_sus.shape)
print("h_res shape:", h_res.shape)
print("pos shape", pos.shape)

# %% compute xpehh
# xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=20000, is_accessible=None, use_threads=True)

# xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, include_edges=True, use_threads=True)
xpehh_raw = allel.xpehh(h_sus, h_res, pos, use_threads=True)
xpehh_raw

# %%  find index of variant with highest xp_ehh signal
xpehh_hit_max = np.nanargmax(xpehh_raw)
xpehh_hit_max

# %% genomic position of top hit
pos[xpehh_hit_max]

# %% Plot raw XP-EHH values

fig, ax = plt.subplots()
ax.hist(xpehh_raw[~np.isnan(xpehh_raw)], bins=20)
ax.set_xlabel('Raw XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize XP-EHH - do not think that we really need to do this

# xpehh_std = allel.standardize_by_allele_count(xpehh_raw, ac_gt[:, 1])
# plot
#fig, ax = plt.subplots()
#ax.hist(xpehh_std[0][~np.isnan(xpehh_std[0])], bins=20)
#ax.set_xlabel('Standardized XP-EHH')
#ax.set_ylabel('Frequency (no. variants)');

# %% look at shapes

print("pos shape:", pos.shape)
print("xpehh_raw shape:", xpehh_raw.shape)

min_pos = pos.min()
max_pos = pos.max()

print("Minimum Genomic Position:", min_pos)
print("Maximum Genomic Position:", max_pos)

# %% plot on manhattan plot

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos, np.abs(xpehh_raw), linestyle=' ', marker='o', mfc='none', mew=3, mec='k', label='$|XP-EHH|$')
ax.axhline(y=4, color='red', linestyle='--')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|XP-EHH|$')
ax.set_ylim(0, 8)
ax.set_xlim(0,61542681)
ax.legend()

# %% list all positions with xpehh value over certain threshold

xpehh_positions_above_threshold_4 = pos[xpehh_raw >= 4]

with open("xpehh_positions_above_threshold_4.txt", "w") as file:
    for position in xpehh_positions_above_threshold_4:
        file.write(str(position) + "\n")

# %% 
## haplotype networks - https://github.com/xgrau/ace1-anopheles-report/blob/master/s01b_haplotype_analysis_Ace1_2020-02-21b.ipynb



# %% H12 was calculated using phased biallelic SNPs in 1000 bp windows
# SNP windows, using the garuds_h function in scikit-allel. 200 permutations?
# Calculate in 1000bp windows, look at the difference in H12
# Calculate for resistant samples

# A.miles:
# To calibrate the window sizes I ran the H12 scans with a range of different window sizes, and chose
# the smallest window size for which the mean value of H1 over all windows was below
# 0.01.
# Lucas et al (2023) to identify regions in which sewpt haplotypes are more frequent in resistant compared to susceptible individuals, they calculated
# the different in H12 value between groups, deltaH12.

h1, h12, h123, h2_h1 = allel.moving_garud_h(h_res_seg, 1000)

# Calculate for susceptible samples

h1, h12, h123, h2_h1 = allel.moving_garud_h(h_sus_seg, 1000)

# look for peaks of difference in H12 between resistant and susceptible samples


# %% PBS
# PBS uses FST to identify genomic regions showing greater evolutionary change in one group 
# (here, the resistant samples) 
# relative to a closely related group (susceptible samples) and an outgroup. While originally designed to detect
# positive selection, it has also been used to detect phenotypic association (Grau-Bové et al., 2021).
# Note, For both H12 and PBS, phenotype permutations were performed as for FST to filter out false positives
# caused by the presence of extended swept haplotypes.
# calculate in 1000 bp windows and plot against the genome. What is classed as significant?
# can I use the control samples as ac3 here?

# %% create ac3 from control samples
# select samples
con_samples = df_samples[df_samples['phenotype'] == 'control'].index.values
con_samples

# select genotypes for samples
gt_con_samples = gt.take(con_samples, axis=1)
gt_con_samples

# create allele counts array
ac_con_samples = gt_con_samples.count_alleles()
ac_con_samples

# %% compute PBS

allel.pbs(ac_sus_samples, ac_res_samples, ac_con_samples, 1000)


### To do:

# make tajima D windows smaller for each chromosome
# work out how to interpret H12 
# look at segregating SNPs for XP-EHH
# output needs to print both chromosome AND position for areas with high iHS and XP-EHH
# linkage disequilibrium?
# hmmIBD (melas)
# GWAS?