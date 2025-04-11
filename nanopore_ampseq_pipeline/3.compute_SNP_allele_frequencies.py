# SNP Frequency Calculator (Filtered by Coverage)
# =====================================================
# This script computes allele frequencies of SNPs from a VCF-derived table,
# taking into account only samples that had sufficient coverage in the
# amplicon where each SNP is located.

import pandas as pd

# ----------------------------
# Configurable file paths
# ----------------------------
snp_file = "/mnt/storage11/sophie/env_dna/double_demux_fastqs/combined_genotyped_filtered_formatted.snps.trans.txt"  # SNPs per sample
tsv_matrix = "amplicon_coverage_matrix.tsv"  # Amplicon x Sample coverage depth matrix
bed_file = "/mnt/storage11/sophie/env_dna/AgamP4_chr.bed"  # BED file with amplicon coordinates
output_file = "snp_frequencies_filtered_by_coverage.tsv"

# ----------------------------
# Step 1: Load input data
# ----------------------------
print("Loading input files...")
snp_df = pd.read_csv(snp_file, sep="\t", header=0, index_col=False)
snp_df["POS"] = snp_df["POS"].astype(int)    
snp_df["CHROM"] = snp_df["CHROM"].astype(str)
coverage_matrix = pd.read_csv(tsv_matrix, sep="\t", index_col=0)
bed_df = pd.read_csv(bed_file, sep="\t", header=None,
                     names=["CHROM", "START", "END", "AMPLICON"])

# ----------------------------
# Step 2: Map SNPs to amplicons
# ----------------------------
def map_snp_to_amplicon(chrom, pos):
    match = bed_df[(bed_df["CHROM"] == chrom) &
                   (bed_df["START"] <= pos) &
                   (bed_df["END"] >= pos)]
    if not match.empty:
        return match.iloc[0]["AMPLICON"]
    return None

print("Mapping SNPs to amplicons...")
snp_df["AMPLICON"] = snp_df.apply(lambda row: map_snp_to_amplicon(row["CHROM"], row["POS"]), axis=1)

# ----------------------------
# Step 3: Convert genotype to ALT allele count
# ----------------------------
def count_alt(gt):
    if pd.isna(gt) or gt == "./.":
        return 0
    if gt == "0/1":
        return 1
    if gt == "1/1":
        return 2
    return 0

snp_df["ALT_COUNT"] = snp_df["GT"].apply(count_alt)

# ----------------------------
# Step 4: Calculate frequencies per SNP
# ----------------------------

print ("Calculating allele frequencies...")
freq_rows = []

for (chrom, pos, alt), group in snp_df.groupby(["CHROM", "POS", "ALT"]):
    amplicon = group["AMPLICON"].iloc[0]
    ref = group["REF"].iloc[0]

    if pd.isna(amplicon) or amplicon not in coverage_matrix.index:
        continue  # Skip if amplicon not found

    # Get samples with sufficient coverage
    covered_samples = coverage_matrix.columns[coverage_matrix.loc[amplicon] > 10]
    subset = group[group["SAMPLE"].isin(covered_samples)]

    alt_sum = subset["ALT_COUNT"].sum()
    n_samples = len(subset)
    freq = alt_sum / (2 * n_samples) if n_samples > 0 else 0

    freq_rows.append({
        "CHROM": chrom,
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "AMPICON": amplicon,
        "ALT_ALLELES": alt_sum,
        "SAMPLES_WITH_DATA": n_samples,
        "FREQ": freq
    })

# ----------------------------
# Step 5: Save output
# ----------------------------
print(f"Writing output to {output_file}...")
freq_df = pd.DataFrame(freq_rows)
freq_df.to_csv(output_file, sep="\t", index=False)

print("Done. SNP frequency table saved.")