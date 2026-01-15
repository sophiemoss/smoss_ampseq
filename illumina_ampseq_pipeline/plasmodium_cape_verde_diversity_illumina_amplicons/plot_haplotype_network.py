#!/usr/bin/env python3
"""
make_haplotype_network.py  —  simplified for replicate-level metadata
(Edited so pie size reflects the number of samples containing the haplotype.)

Build a haplotype network (MST) with pie-chart nodes from:
  1) FASTA of unique, aligned haplotype sequences (equal length), and
  2) a counts table (HaplotypR-style) with per-sample read counts,
  3) optional metadata joined directly on SampleID (e.g., Year).

Outputs
-------
- <OUT>_network.png    : PNG with network + legend embedded
- <OUT>_network.gexf   : graph file for Gephi/Cytoscape
- <OUT>_node_table.csv : node totals + per-group fractions
- <OUT>_edge_table.csv : (opt) edge table with Hamming distances

Install (once):
  pip install biopython pandas numpy matplotlib networkx

Example:
python /mnt/storage11/sophie/gitrepos/smoss_ampseq/illumina_ampseq_pipeline/plasmodium_cape_verde_diversity_illumina_amplicons/plot_haplotype_network.py \
  --fasta "ama1_HaplotypeSeq_merge.fasta" \
  --counts "AMA1_finalTab_passed_technical_replicates.txt" \
  --metadata "AMA1_samplespassedtechnicalreps.csv" \
  --group Year \
  --marker AMA1 \
  --out-prefix ama1_network_by_year \
  --save-edges

Notes
-----
- Node size is now the number of unique SampleIDs containing each haplotype.
- Pie slices still reflect fractions by the --group column (e.g., Year).
"""

import argparse
from pathlib import Path
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from Bio import SeqIO
from matplotlib.patches import Wedge, Patch

ARTEFACTS = {"Noise", "Indels", "Singelton", "Chimera", "NA"}  # labels to drop


# --------------------- Robust readers / headers ---------------------

def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Strip BOM/whitespace; alias common names case-insensitively."""
    cols = [str(c).replace("\ufeff", "").strip() for c in df.columns]
    df.columns = cols
    lc = {c.lower(): c for c in df.columns}
    if "sampleid" in lc and "SampleID" not in df.columns:
        df.rename(columns={lc["sampleid"]: "SampleID"}, inplace=True)
    if "haplotype" in lc and "Haplotype" not in df.columns:
        df.rename(columns={lc["haplotype"]: "Haplotype"}, inplace=True)
    if "reads" in lc and "Reads" not in df.columns:
        df.rename(columns={lc["reads"]: "Reads"}, inplace=True)
    return df

def _read_table_auto(path: str) -> pd.DataFrame:
    """Auto-detect delimiter, then normalize headers."""
    df = pd.read_csv(path, sep=None, engine="python")
    return _normalize_columns(df)


# ------------------------- I/O helpers -------------------------

def read_haplotypes_fasta(fasta_path: str) -> Dict[str, str]:
    """Load haplotype sequences; uppercase, strip gaps/spaces; enforce equal length."""
    seqs: Dict[str, str] = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        sid = str(rec.id).strip()
        seq = str(rec.seq).upper().replace("-", "").replace(" ", "")
        seqs[sid] = seq
    if not seqs:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")
    lens = {len(s) for s in seqs.values()}
    if len(lens) != 1:
        raise ValueError(f"Sequences are not equal length (lengths: {sorted(lens)}). Provide aligned haplotypes.")
    return seqs


def load_counts_table(counts_path: str) -> pd.DataFrame:
    """Load counts (SampleID, Haplotype, Reads); filter artefacts and zeros."""
    df = _read_table_auto(counts_path)
    required = {"SampleID", "Haplotype", "Reads"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Counts file missing columns: {missing} (found: {list(df.columns)})")
    df = df.copy()
    df = df[~df["Haplotype"].isin(ARTEFACTS)]
    df["Reads"] = pd.to_numeric(df["Reads"], errors="coerce").fillna(0).astype(float)
    df = df[df["Reads"] > 0]
    return df[["SampleID", "Haplotype", "Reads"]]


def merge_metadata_simple(counts: pd.DataFrame, metadata_path: str) -> pd.DataFrame:
    """
    Simple merge of metadata on SampleID (your metadata is replicate-level).
    The metadata file must contain 'SampleID' and the grouping column(s).
    """
    meta = _read_table_auto(metadata_path)
    if "SampleID" not in meta.columns:
        raise ValueError(f"Metadata must contain a 'SampleID' column. Found: {list(meta.columns)}")
    return counts.merge(meta, on="SampleID", how="left")


# ------------------------- Graph construction -------------------------

def hamming(a: str, b: str) -> int:
    """Simple Hamming distance between equal-length strings."""
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def build_mst(seqs: Dict[str, str]) -> nx.Graph:
    """Complete graph with Hamming weights → minimum spanning tree (Kruskal)."""
    ids = list(seqs.keys())
    G = nx.Graph()
    G.add_nodes_from(ids)
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            d = hamming(seqs[ids[i]], seqs[ids[j]])
            G.add_edge(ids[i], ids[j], weight=d)
    return nx.minimum_spanning_tree(G, weight="weight", algorithm="kruskal")


def pie_fracs_for_nodes(counts: pd.DataFrame, hap_ids: List[str], group_col: str, presence: bool = False
                        ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return (counts_wide, fracs_wide): rows=Haplotype, columns=groups.
    presence=True → presence/absence per group; else sum of Reads per group.
    """
    if group_col not in counts.columns:
        raise ValueError(f"Grouping column '{group_col}' not found. Available: {sorted(counts.columns)}")

    if presence:
        agg = (counts.assign(present=1)
                      .groupby(["Haplotype", group_col])["present"]
                      .max()
                      .reset_index()
                      .rename(columns={"present": "value"}))
    else:
        agg = (counts.groupby(["Haplotype", group_col])["Reads"]
                     .sum()
                     .reset_index()
                     .rename(columns={"Reads": "value"}))

    counts_wide = agg.pivot(index="Haplotype", columns=group_col, values="value").fillna(0.0)
    counts_wide = counts_wide.reindex(hap_ids).fillna(0.0)
    totals = counts_wide.sum(axis=1).replace(0, np.nan)
    fracs_wide = counts_wide.div(totals, axis=0).fillna(0.0)
    return counts_wide, fracs_wide


# ------------------------- Drawing -------------------------

def build_color_list(n: int) -> List[str]:
    """Matplotlib default cycle (cycled as needed)."""
    base = plt.rcParams['axes.prop_cycle'].by_key().get('color', [])
    if not base:
        base = ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"]
    return [base[i % len(base)] for i in range(n)]

def draw_pie(ax, center, r, fracs, colors: List[str]):
    """Draw a pie at center with radius r using fracs + colors."""
    start = 0.0
    for frac, color in zip(fracs, colors):
        if frac <= 0:
            continue
        theta = 360.0 * float(frac)
        wedge = Wedge(center, r, start, start + theta, linewidth=0)
        wedge.set_facecolor(color)
        ax.add_patch(wedge)
        start += theta


def plot_network(MST: nx.Graph,
                 counts_wide: pd.DataFrame,
                 fracs_wide: pd.DataFrame,
                 size_totals: pd.Series,     # NEW: node-size basis = #unique SampleIDs per haplotype
                 out_png: str,
                 title: str = "",
                 layout: str = "spring",
                 min_r: float = 0.03,
                 max_r: float = 0.10,
                 legend_loc: str = "upper right",
                 legend_ncol: int = 1):
    """Plot MST with pie nodes; legend embedded on same figure."""
    # Layout
    if layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(MST)
    elif layout == "circular":
        pos = nx.circular_layout(MST)
    else:
        pos = nx.spring_layout(MST, seed=42)

    # === Colors (by group names) ===
    group_names = list(fracs_wide.columns)
    custom_colors = {
        "2024": "mediumslateblue",
        "2025": "darkorange",
    }
    base_palette = build_color_list(len(group_names))
    colors = [custom_colors.get(str(g), base_palette[i]) for i, g in enumerate(group_names)]

    fig, ax = plt.subplots(figsize=(10, 8))

    # Edges
    for u, v, d in MST.edges(data=True):
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], '-', lw=1, color="black", alpha=0.6)

    # --- Node radii from size_totals (number of unique SampleIDs per haplotype) ---
    # Ensure index alignment to haplotypes in fracs_wide
    totals = size_totals.reindex(fracs_wide.index).fillna(0.0)
    if (totals.max() - totals.min()) <= 0:
        radii = pd.Series(min_r, index=totals.index)
    else:
        radii = min_r + (totals - totals.min()) / (totals.max() - totals.min()) * (max_r - min_r)

    # Nodes (pies + outlines + labels)
    for hap in MST.nodes():
        x, y = pos[hap]
        r = radii.loc[hap] if hap in radii.index else min_r
        fracs = fracs_wide.loc[hap].values if hap in fracs_wide.index else np.array([1.0])
        s = fracs.sum()
        fracs = (fracs / s) if s > 0 else np.array([1.0])
        draw_pie(ax, (x, y), r, fracs, colors)
        ax.add_patch(plt.Circle((x, y), r, fill=False, linewidth=0.9, color="black"))
        ax.text(x, y + r + 0.02, hap, ha='center', va='bottom', fontsize=8)

    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title or "Haplotype network (MST)")

    # Embedded legend
    if len(group_names) > 0:
        proxies = [Patch(facecolor=c, edgecolor="none") for c in colors]
        labels  = [str(n) for n in group_names]
        leg = ax.legend(proxies, labels, loc=legend_loc, ncol=legend_ncol,
                        frameon=True, framealpha=0.9, title="Group")
        leg._legend_box.align = "left"

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------- CLI -------------------------

def save_gexf(G: nx.Graph, out_gexf: str):
    nx.write_gexf(G, out_gexf)

def main():
    ap = argparse.ArgumentParser(
        description="Build a haplotype MST network (pie nodes) from FASTA + counts, with optional metadata (e.g., Year)."
    )
    ap.add_argument("--fasta", required=True, help="Haplotype FASTA (e.g., csp_HaplotypeSeq_merge.fasta)")
    ap.add_argument("--counts", required=True, help="Counts table (TSV/CSV) with at least SampleID,Haplotype,Reads")
    ap.add_argument("--metadata", default=None, help="Optional metadata CSV/TSV with SampleID and extra fields (e.g., Year, Site).")
    ap.add_argument("--group", default="SampleID",
                    help="Grouping column to color pies by (e.g., Year/Site/SampleID).")
    ap.add_argument("--marker", default="", help="Marker label for figure title (e.g., CSP or TRAP)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix for files")
    ap.add_argument("--layout", choices=["spring","kamada_kawai","circular"], default="spring",
                    help="Graph layout algorithm for plotting (default: spring)")
    ap.add_argument("--presence", action="store_true",
                    help="Use presence/absence per group for pie *fractions*. (Node size is by #samples.)")
    ap.add_argument("--min-radius", type=float, default=0.03, help="Minimum node radius (plot units).")
    ap.add_argument("--max-radius", type=float, default=0.10, help="Maximum node radius (plot units).")
    ap.add_argument("--legend-loc", default="upper right", help="Legend location (matplotlib string).")
    ap.add_argument("--legend-ncol", type=int, default=1, help="Legend columns.")
    ap.add_argument("--save-edges", action="store_true", help="Also save an edge table CSV with Hamming distances.")
    args = ap.parse_args()

    # Read inputs
    seqs = read_haplotypes_fasta(args.fasta)
    counts = load_counts_table(args.counts)

    # Merge metadata directly on SampleID (if provided)
    if args.metadata:
        counts = merge_metadata_simple(counts, args.metadata)

    # Ensure requested group exists and has values
    if args.group not in counts.columns:
        raise ValueError(f"Requested group '{args.group}' not found after merge. "
                         f"Available columns: {sorted(counts.columns)}")
    if counts[args.group].dropna().empty:
        print(f"[WARN] Column '{args.group}' has no values after merge; legend/pies will be uncolored.")

    # Keep only haplotypes present in the FASTA
    counts = counts[counts["Haplotype"].isin(seqs.keys())].copy()

    # Keep only haplotypes present in the counts (accepted haplotypes that passed filtering)
    seqs = {k: v for k, v in seqs.items() if k in counts["Haplotype"].unique()}

    # Debug lists
    print("Haplotypes in counts:", sorted(counts["Haplotype"].unique()))
    print("Haplotypes in filtered FASTA:", sorted(seqs.keys()))

    # Build MST
    MST = build_mst(seqs)

    # Compute pies (fractions by group)
    hap_ids = list(MST.nodes())
    counts_wide, fracs_wide = pie_fracs_for_nodes(counts, hap_ids, group_col=args.group, presence=args.presence)

    # --- Compute node-size totals: #unique SampleIDs per haplotype ---
    size_totals = (counts.drop_duplicates(["Haplotype", "SampleID"])
                          .groupby("Haplotype")["SampleID"]
                          .nunique())

    # Outputs
    out_png   = f"{args.out_prefix}_network.png"
    out_gexf  = f"{args.out_prefix}_network.gexf"
    out_nodes = f"{args.out_prefix}_node_table.csv"
    out_edges = f"{args.out_prefix}_edge_table.csv"

    # Plot
    plot_network(MST, counts_wide, fracs_wide, size_totals, out_png,
                 title=f"{args.marker} Haplotype Network (MST)" if args.marker else "Haplotype Network (MST)",
                 layout=args.layout,
                 min_r=args.min_radius,
                 max_r=args.max_radius,
                 legend_loc=args.legend_loc,
                 legend_ncol=args.legend_ncol)

    # Save graph + node table (total_value = #samples per haplotype)
    totals = size_totals.reindex(fracs_wide.index).fillna(0.0).rename("total_value")
    node_table = pd.concat([totals, fracs_wide], axis=1)
    node_table.to_csv(out_nodes, index=True)

    # Optional: edge table
    if args.save_edges:
        rows = [{"source": u, "target": v, "hamming": d.get("weight", np.nan)} for u, v, d in MST.edges(data=True)]
        pd.DataFrame(rows).to_csv(out_edges, index=False)

    print(f"Wrote: {out_png}")
    print(f"Wrote: {out_gexf}")
    print(f"Wrote: {out_nodes}")
    if args.save_edges:
        print(f"Wrote: {out_edges}")


if __name__ == "__main__":
    main()
