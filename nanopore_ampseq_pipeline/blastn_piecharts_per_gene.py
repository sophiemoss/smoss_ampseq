import pandas as pd
import matplotlib.pyplot as plt
import re
from collections import defaultdict, Counter
import glob
import os

# === Custom color assignment for species ===
species_colors = {
    'Anopheles gambiae': '#59c3c3',
    'Anopheles coluzzii': '#52489c',
    # Add more species and colors here
}
default_color = '#bbbbbb'  # fallback for species not explicitly colored


# === Get input files ===
files = glob.glob("sample*_species_markers_merged_nt_results.txt")

# === Count species per gene across all files ===
gene_species_counts = defaultdict(Counter)
gene_sample_counts = defaultdict(set)  # track samples with ≥5 hits per gene

for file in files:
    sample_name = os.path.basename(file).split("_")[0]  # e.g., "sample10"
    gene_hits = defaultdict(list)

    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            query = parts[0].strip()
            stitle = parts[-1].strip()

            gene_match = re.search(r'_(\w+)$', query)
            if not gene_match:
                continue
            gene = gene_match.group(1)

            species_match = re.search(r'([A-Z][a-z]+ [a-z0-9\-]+)', stitle)
            if not species_match:
                continue
            species = species_match.group(1)

            gene_species_counts[gene][species] += 1
            gene_hits[gene].append(species)

    # After processing the file, check how many genes got ≥5 hits
    for gene, hits in gene_hits.items():
        if len(hits) >= 5:
            gene_sample_counts[gene].add(sample_name)

# === Plot pie chart per gene ===
for gene, species_counter in gene_species_counts.items():
    df = pd.DataFrame(species_counter.items(), columns=["Species", "Hits"])
    df["%"] = (df["Hits"] / df["Hits"].sum()) * 100

    # Color setup
    pie_colors = [species_colors.get(species, default_color) for species in df["Species"]]
    sample_count = len(gene_sample_counts[gene])

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, texts, autotexts = ax.pie(
        df["%"], labels=None, colors=pie_colors,
        autopct="%1.1f%%", startangle=90, textprops={'fontsize': 10}
    )
    ax.set_title(f"{gene} Blastn Species Composition", fontsize=14)

    # Create custom legend entries
    legend_labels = [f"{species}" for species in df["Species"]]
    legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=species_colors.get(species, default_color), markersize=10)
                      for species in df["Species"]]

    # Add 'n = X' as the first legend entry
    legend_labels.insert(0, f"N = {sample_count} samples")
    legend_handles.insert(0, plt.Line2D([0], [0], color='none'))  # invisible entry

    # Add legend outside plot
    ax.legend(
        legend_handles,
        legend_labels,
        title="Species",
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        fontsize=9,
        title_fontsize=10
    )

    plt.tight_layout()
    filename = f"{gene}_species_piechart.png"
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.close()

    print(f"✅ Saved pie chart for {gene} with legend (n = {sample_count}): {filename}")
