import pandas as pd
import matplotlib.pyplot as plt
import re
from collections import Counter, defaultdict
import glob
import os
from datetime import datetime

# === Logging setup ===
log_file = "blastn_species_summary_log.txt"
def log(message):
    print(message)
    with open(log_file, 'a') as f:
        f.write(message + '\n')

# === Start new log ===
with open(log_file, 'w') as f:
    f.write(f"📅 Run started: {datetime.now()}\n")

# === Custom color assignment for species ===
species_colors = {
    'Anopheles gambiae': '#59c3c3',
    'Anopheles coluzzii': '#52489c',
    # Add more species and colors here
}
default_color = '#bbbbbb'  # fallback for species not explicitly colored

# === Load all *_nt_results.txt files ===
result_files = glob.glob("sample*_species_markers_merged_nt_results.txt")
log(f"🗂 Found {len(result_files)} result file(s).\n")

# === Parse files ===
data = defaultdict(Counter)

for file_path in result_files:
    log(f"🔍 Processing: {file_path}")
    species_hits = []
    empty_lines = 0
    malformed_lines = 0
    unmatched_lines = 0

    with open(file_path, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                empty_lines += 1
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                malformed_lines += 1
                continue

            stitle = parts[-1].strip()
            match = re.search(r'([A-Z][a-z]+ [a-z0-9\-]+)', stitle)
            if match:
                species = match.group(1)
                species_hits.append(species)
            else:
                unmatched_lines += 1
                log(f"⚠️ No species match at line {line_number}: {stitle}")

    if not species_hits:
        log(f"🚫 No valid species hits found in {file_path}, skipping.\n")
        continue

    sample_name_match = re.match(r'(sample\d+)_species_markers_merged_nt_results\.txt', os.path.basename(file_path))
    sample_name = sample_name_match.group(1) if sample_name_match else os.path.basename(file_path)

    total_hits = len(species_hits)
    counts = Counter(species_hits)
    for species, count in counts.items():
        data[species][sample_name] = (count / total_hits) * 100  # % of total hits

    log(f"✅ Parsed {sample_name}: {total_hits} hits, {empty_lines} empty, {malformed_lines} malformed, {unmatched_lines} unmatched\n")

# === Build dataframe: species x sample ===
df = pd.DataFrame(data).fillna(0).T
df = df.sort_index(axis=1)  # sort samples

# === Plot stacked bar chart with manual species colors ===
fig, ax = plt.subplots(figsize=(12, 8))
bottom = [0] * len(df.columns)
species_list = list(df.index)

for species in species_list:
    color = species_colors.get(species, default_color)
    # italicise the species names while retaining the space between genus and species!
    genus_species = species.replace(" ", r"\ ")
    italic_species = f"$\\it{{{genus_species}}}$"
    ax.bar(df.columns, df.loc[species], bottom=bottom, label=italic_species, color=color)
    bottom = [i + j for i, j in zip(bottom, df.loc[species])]

ax.set_ylabel('% of Total Hits')
ax.set_xlabel('Sample')
ax.set_title('Species Composition by Sample (BLAST top 5 hits per gene)')
ax.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("blastn_species_stackedbarchart.png", dpi=600)
plt.close()

log("📊 Stacked bar chart saved as: blastn_species_stackedbarchart.png")
log("🏁 All processing complete.\n")
