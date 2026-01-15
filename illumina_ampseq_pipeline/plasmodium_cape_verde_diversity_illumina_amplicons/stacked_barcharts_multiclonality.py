#!/usr/bin/env python3
# stacked_bars_multiclonality.py
#
# Build stacked bar charts for multiclonality:
#  - By YEAR (reads + proportions)
#  - By SAMPLE or SPECIMEN (reads + proportions, with year labels on bars)
#
"""
Usage example:

python /mnt/storage11/sophie/gitrepos/smoss_ampseq/illumina_ampseq_pipeline/plasmodium_cape_verde_diversity_illumina_amplicons/stacked_barcharts_multiclonality.py \
    --infile CSP_finalTab_passed_technical_replicates.txt \
    --marker csp \
    --outdir new_figures \
    --sep "\t" \
    --year-map CSP_samplespassedtechnicalreps.csv \
    --year-map-sep "," \
    --dpi 600 \
    --collapse-replicates
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from statistics import mode


def parse_args():
    p = argparse.ArgumentParser(
        description="Stacked bars of haplotypes (by year and by sample/specimen) for a given marker."
    )
    p.add_argument("--infile", required=True, help="Input TSV/CSV with haplotype rows.")
    p.add_argument("--sep", default="\t", help=r"Field separator for infile (default: '\t').")
    p.add_argument("--marker", required=True, help="Marker to filter (e.g., csp, trap, ama1).")
    p.add_argument("--outdir", default="figures", help="Output directory (default: figures).")
    p.add_argument("--year-column", default="Year", help="Name of year column if it exists (default: Year).")
    p.add_argument("--year-map", help="Optional mapping (SampleID,Year) if Year column is missing.")
    p.add_argument("--year-map-sep", default=",", help=r"Separator for year-map (default ','). Use '\t' for tab.")
    p.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Resolution for output figures (default: 600 dpi, journal quality)"
    )
    p.add_argument(
        "--collapse-replicates",
        action="store_true",
        help="Collapse technical replicates to one bar per SPECIMEN in per-sample plots."
    )
    return p.parse_args()


def unescape_sep(s: str) -> str:
    """Turn literal sequences like '\\t' or '\\n' into actual tab/newline."""
    try:
        return s.encode("utf-8").decode("unicode_escape")
    except Exception:
        return s


def ensure_year(df, year_col, year_map_path, year_map_sep):
    if year_col in df.columns:
        out = df.rename(columns={year_col: "Year"})
    elif year_map_path:
        ym = pd.read_csv(year_map_path, sep=year_map_sep)
        if not {"SampleID", "Year"}.issubset(ym.columns):
            raise ValueError("Year map must contain columns: SampleID, Year")
        out = df.merge(ym[["SampleID", "Year"]], on="SampleID", how="left")
    else:
        raise RuntimeError("Need a Year column in the infile or provide --year-map")

    if out["Year"].isna().any():
        missing = out.loc[out["Year"].isna(), "SampleID"].unique()
        raise ValueError(
            f"Year is missing for {len(missing)} SampleID(s). "
            f"Examples: {', '.join(map(str, missing[:5]))} ..."
        )
    return out


def stacked_bar(df_wide, x_labels, title, outfile, args,
                annotate_years=None, ylabel="Total reads",
                is_sample_plot=False, x_label_override=None):

    # Taller for sample/specimen plots
    height = 10 if is_sample_plot else 8
    fig, ax = plt.subplots(figsize=(max(10, 0.5 * len(x_labels)), height))

    bottom = None
    for hap in sorted(df_wide.columns):
        ax.bar(x_labels, df_wide[hap], bottom=bottom, label=hap)
        bottom = df_wide[hap] if bottom is None else (bottom + df_wide[hap])

    # --- unchanged: annotate reads or proportions ---
    if annotate_years is not None:
        if "proportion" in ylabel.lower():
            # proportions: annotate above bars
            for ix, lab in enumerate(x_labels):
                ax.text(ix, 1.02, str(annotate_years.loc[lab]),
                        rotation=90, va="bottom", ha="center", fontsize=7)
        else:
            totals = df_wide.sum(axis=1)
            for ix, (lab, total) in enumerate(zip(x_labels, totals)):
                y = total + max(1, total * 0.01)
                ax.text(ix, y, str(annotate_years.loc[lab]),
                        rotation=90, va="bottom", ha="center", fontsize=7)

    # Labels
    ax.set_ylabel(ylabel)
    xlab = x_label_override if x_label_override else ("Year" if not is_sample_plot else "SampleID")
    ax.set_xlabel(xlab)

    # ✅ Put title outside plotting area
    fig.suptitle(title, y=0.98, fontsize=14)   # adjust y if needed

    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=90 if is_sample_plot else 0)

    ax.legend(title="Haplotype", bbox_to_anchor=(1.02, 1), loc="upper left")

    # ✅ layout: first tighten, then reserve top margin
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)   # leaves space for suptitle

    fig.savefig(outfile, dpi=args.dpi)
    plt.close(fig)
    print(f"[✔] Saved: {outfile}")


def main():
    args = parse_args()
    args.sep = unescape_sep(args.sep)
    args.year_map_sep = unescape_sep(args.year_map_sep)

    marker = args.marker.lower()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Plotting {marker.upper()} haplotypes ===")
    print(f"Input table : {args.infile}")
    if args.year_map:
        print(f"Year map    : {args.year_map} (sep='{args.year_map_sep}')")
    print(f"Output dir  : {outdir}")
    print(f"DPI         : {args.dpi}")
    print(f"Collapse reps: {'YES' if args.collapse_replicates else 'NO'}\n")

    # Load input and keep only what we need
    df = pd.read_csv(args.infile, sep=args.sep)
    needed = ["SampleID", "SampleName", "MarkerID", "Haplotype", "Reads", "Specimen"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Infile missing columns: {missing}")

    # Subset and filter to marker
    df = df[needed + ([args.year_column] if args.year_column in df.columns else [])].copy()
    df = df[df["MarkerID"].astype(str).str.lower() == marker]

    # Ensure Year column exists/merged
    df = ensure_year(df, args.year_column, args.year_map, args.year_map_sep)

    # -----------------------------
    # STACKED BY YEAR (reads)
    # -----------------------------
    by_year = df.groupby(["Year", "Haplotype"], as_index=False)["Reads"].sum()
    pivot_year = by_year.pivot(index="Year", columns="Haplotype", values="Reads").fillna(0).sort_index()

    stacked_bar(
        df_wide=pivot_year,
        x_labels=pivot_year.index.astype(str),
        title=f"{marker.upper()} haplotypes per year (stacked reads)",
        outfile=outdir / f"{marker}_stacked_by_year_reads.png",
        args=args,
        is_sample_plot=False,
        x_label_override="Year"
    )

    # Proportions by year
    prop_year = pivot_year.div(pivot_year.sum(axis=1).replace(0, 1), axis=0)
    stacked_bar(
        df_wide=prop_year,
        x_labels=prop_year.index.astype(str),
        title=f"{marker.upper()} haplotypes per year (stacked proportions)",
        outfile=outdir / f"{marker}_stacked_by_year_proportions.png",
        args=args,
        ylabel="Proportion",
        is_sample_plot=False,
        x_label_override="Year"
    )

    # ---------------------------------------------------------
    # STACKED BY SAMPLE or SPECIMEN - COLLAPSING THE REPLICATES
    # ---------------------------------------------------------
    if args.collapse_replicates:
        # Collapse to SPECIMEN: sum reads across SampleIDs within each Specimen + haplotype
        # Derive a Year per SPECIMEN as mode across its samples (ties → min)
        sid_year = df[["SampleID", "Year"]].drop_duplicates()
        spec_map = df[["SampleID", "Specimen"]].drop_duplicates()
        spec_years = sid_year.merge(spec_map, on="SampleID").groupby("Specimen")["Year"].apply(list)

        def mode_or_min(lst):
            try:
                return mode(lst)
            except Exception:
                return min(lst)

        year_per_specimen = spec_years.apply(mode_or_min)  # index = Specimen

        # Merge specimen-level year into df; avoid Year_x/Year_y by using a new column
        df_spec = df.merge(year_per_specimen.rename("Year_specimen"),
                           left_on="Specimen", right_index=True, how="left")

        # Choose the plotting year: prefer specimen-level; fall back to row Year
        if "Year" in df_spec.columns:
            df_spec["Year_plot"] = df_spec["Year_specimen"].fillna(df_spec["Year"])
        else:
            df_spec["Year_plot"] = df_spec["Year_specimen"]

        if df_spec["Year_plot"].isna().any():
            missing_specs = df_spec.loc[df_spec["Year_plot"].isna(), "Specimen"].unique()
            print(f"[WARN] Missing Year for {len(missing_specs)} specimen(s): {missing_specs[:5]} ...")

        # Aggregate reads at Specimen × Year_plot × Haplotype
        by_specimen = (df_spec
                       .groupby(["Specimen", "Year_plot", "Haplotype"], as_index=False)["Reads"]
                       .sum())

        pivot_sample = (by_specimen
                        .pivot(index="Specimen", columns="Haplotype", values="Reads")
                        .fillna(0))

        # Year labels aligned to plotted index
        year_axis_series = (by_specimen.drop_duplicates(subset=["Specimen"])
                            .set_index("Specimen")["Year_plot"]
                            .reindex(pivot_sample.index))

        # Outputs and labels
        sample_reads_png = outdir / f"{marker}_stacked_by_specimen_reads_year_labeled.png"
        sample_prop_png  = outdir / f"{marker}_stacked_by_specimen_proportions_year_labeled.png"
        sample_title_reads = f"{marker.upper()} haplotypes per specimen (stacked reads); year annotated"
        sample_title_prop  = f"{marker.upper()} haplotypes per specimen (stacked proportions); year annotated"
        x_label = "Specimen"
        is_sample_plot = True

    else:
        # Keep SampleID (Rep1/Rep2 visible)
        by_sample = df.groupby(["SampleID", "Year", "Haplotype"], as_index=False)["Reads"].sum()
        pivot_sample = by_sample.pivot(index="SampleID", columns="Haplotype", values="Reads").fillna(0)
        year_axis_series = by_sample.drop_duplicates(subset=["SampleID"]).set_index("SampleID")["Year"]

        sample_reads_png = outdir / f"{marker}_stacked_by_sample_reads_year_labeled.png"
        sample_prop_png  = outdir / f"{marker}_stacked_by_sample_proportions_year_labeled.png"
        sample_title_reads = f"{marker.upper()} haplotypes per sample (stacked reads); year annotated"
        sample_title_prop  = f"{marker.upper()} haplotypes per sample (stacked proportions); year annotated"
        x_label = "SampleID"
        is_sample_plot = True

    # Reads plot (sample/specimen)
    stacked_bar(
        df_wide=pivot_sample,
        x_labels=pivot_sample.index,
        title=sample_title_reads,
        outfile=sample_reads_png,
        args=args,
        annotate_years=year_axis_series,
        ylabel="Reads",
        is_sample_plot=is_sample_plot,
        x_label_override=x_label
    )

    # Proportions plot (sample/specimen)
    prop_sample = pivot_sample.div(pivot_sample.sum(axis=1).replace(0, 1), axis=0)
    stacked_bar(
        df_wide=prop_sample,
        x_labels=prop_sample.index,
        title=sample_title_prop,
        outfile=sample_prop_png,
        args=args,
        annotate_years=year_axis_series,
        ylabel="Proportion",
        is_sample_plot=is_sample_plot,
        x_label_override=x_label
    )

    print("\n✅ All plots complete.")
    print(f"Outputs saved in: {outdir}\n")


if __name__ == "__main__":
    main()
