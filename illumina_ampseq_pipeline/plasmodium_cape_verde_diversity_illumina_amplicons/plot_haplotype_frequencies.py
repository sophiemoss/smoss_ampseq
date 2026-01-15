#!/usr/bin/env python3
# haplotype_frequency.py
#
# Plot haplotype frequencies for a chosen marker:
#   1) # of samples carrying each haplotype
#   2) Total reads per haplotype
#

'''
python /mnt/storage11/sophie/gitrepos/smoss_ampseq/illumina_ampseq_pipeline/plasmodium_cape_verde_diversity_illumina_amplicons/plot_haplotype_frequencies.py \
  --infile TRAP_finalTab_passed_technical_replicates.txt \
  --marker trap \
  --outdir new_figures \
  --sep "\t" \
  --dpi 600
'''


import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Plot haplotype frequencies for a chosen marker.")
    p.add_argument("--infile", required=True, help="Input TSV/CSV haplotype table.")
    p.add_argument("--sep", default="\t", help=r"Field separator for infile (default: '\t').")
    p.add_argument("--marker", required=True, help="Marker to filter (e.g., csp, trap, ama1).")
    p.add_argument("--marker-column", default="MarkerID", help="Column name that holds the marker (default: MarkerID).")
    p.add_argument("--hap-column", default="Haplotype", help="Column name for haplotype (default: Haplotype).")
    p.add_argument("--reads-column", default="Reads", help="Column name for reads (default: Reads).")
    p.add_argument("--sampleid-column", default="SampleID", help="Column name for sample ID (default: SampleID).")
    p.add_argument("--outdir", default="figures", help="Output directory (default: figures).")
    p.add_argument("--dpi", type=int, default=600, help="Figure resolution (default: 600).")
    return p.parse_args()


def unescape_sep(s: str) -> str:
    try:
        return s.encode("utf-8").decode("unicode_escape")
    except Exception:
        return s


def barplot(series: pd.Series, title: str, outfile: Path, dpi: int):
    # Series index = haplotypes, values = counts
    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(series.index)), 6))

    ax.bar(series.index.astype(str), series.values)
    ax.set_ylabel("Count")
    ax.set_xlabel("Haplotype")
    ax.set_xticklabels(series.index.astype(str), rotation=45, ha="right")

    # Put title outside axes to avoid overlap
    fig.suptitle(title, y=0.98, fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)
    print(f"[✔] Saved: {outfile}")


def main():
    args = parse_args()
    args.sep = unescape_sep(args.sep)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Haplotype frequencies for {args.marker.upper()} ===")
    print(f"Input  : {args.infile}")
    print(f"Outdir : {outdir}")
    print(f"DPI    : {args.dpi}\n")

    # Load & basic checks
    df = pd.read_csv(args.infile, sep=args.sep)

    for col in [args.marker_column, args.hap_column, args.reads_column, args.sampleid_column]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in input table.")

    # Filter marker
    marker = str(args.marker).lower()
    df = df[df[args.marker_column].astype(str).str.lower() == marker].copy()

    if df.empty:
        raise SystemExit(f"No rows found for marker '{args.marker}' in {args.infile}")

    # Compute frequencies
    # 1) Number of samples carrying each haplotype
    hap_samples = (
        df.groupby(args.hap_column)[args.sampleid_column]
          .nunique()
          .sort_values(ascending=False)
    )

    # 2) Total reads per haplotype
    hap_reads = (
        df.groupby(args.hap_column)[args.reads_column]
          .sum()
          .sort_values(ascending=False)
    )

    # Plots
    out1 = outdir / f"{marker}_haplotype_frequency_samples.png"
    barplot(
        hap_samples,
        title=f"{args.marker.upper()} haplotype frequency (samples)",
        outfile=out1,
        dpi=args.dpi,
    )

    out2 = outdir / f"{marker}_haplotype_frequency_reads.png"
    # Reuse barplot but adjust y label after
    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(hap_reads.index)), 6))
    ax.bar(hap_reads.index.astype(str), hap_reads.values)
    ax.set_ylabel("Total reads")
    ax.set_xlabel("Haplotype")
    ax.set_xticklabels(hap_reads.index.astype(str), rotation=45, ha="right")
    fig.suptitle(f"{args.marker.upper()} haplotype frequency (reads)", y=0.98, fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    fig.savefig(out2, dpi=args.dpi)
    plt.close(fig)
    print(f"[✔] Saved: {out2}")

    # Optional: export a small table
    freq_table = pd.DataFrame({
        "samples": hap_samples,
        "reads": hap_reads.reindex(hap_samples.index).fillna(0).astype(int)
    })
    out3 = outdir / f"{marker}_haplotype_frequency_table.tsv"
    freq_table.to_csv(out3, sep="\t")
    print(f"[✔] Saved: {out3}\n")

    print("✅ Done.")


if __name__ == "__main__":
    main()
