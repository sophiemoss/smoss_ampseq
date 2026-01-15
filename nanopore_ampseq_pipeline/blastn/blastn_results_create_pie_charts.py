import pandas as pd
import matplotlib.pyplot as plt
import re
from collections import Counter
import glob
import os
from datetime import datetime

# === Setup log file ===
log_file = "blastn_species_summary_log.txt"
log_mode = 'w'  # Overwrite each time. Use 'a' to append instead.

def log(message):
    print(message)
    with open(log_file, log_mode if not os.path.exists(log_file) else 'a') as f:
        f.write(message + '\n')

# === Step 1: Find all *_nt_results.txt files ===
result_files = glob.glob("sample*_species_markers_merged_nt_results.txt")
log(f"\n📅 Run started: {datetime.now()}\n")
log(f"🗂 Found {len(result_files)} result file(s) to process.")

# === Step 2: Process each file ===
for file_path in result_files:
    log(f"\n🔍 Processing: {file_path}")
    
    species_hits = []
    empty_lines = 0
    malformed_lines = 0
    unmatched_lines = 0

    with open(file_path, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                empty_lines += 1
                log(f"⚠️ Empty line at line {line_number} in {file_path}")
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                malformed_lines += 1
                log(f"⚠️ Malformed line (only {len(parts)} fields) at line {line_number} in {file_path}")
                continue

            # Parse the final column for species info
            stitle = parts[-1].strip()

            # Extract species name (Genus species), allowing lowercase species + optional hyphen
            match = re.search(r'([A-Z][a-z]+ [a-z0-9\-]+)', stitle)
            if match:
                species = match.group(1)
                species_hits.append(species)
            else:
                unmatched_lines += 1
                log(f"⚠️ No species match in stitle at line {line_number}: {stitle}")

    if not species_hits:
        log(f"🚫 No valid species hits found in {file_path}, skipping chart.")
        continue

    # Create summary DataFrame
    species_counts = Counter(species_hits)
    df = pd.DataFrame(species_counts.items(), columns=['Species', 'Hit Count'])
    df['% of Total Hits'] = (df['Hit Count'] / len(species_hits)) * 100

    # Extract sample name
    sample_name_match = re.match(r'(sample\d+)_species_markers_merged_nt_results\.txt', os.path.basename(file_path))
    sample_name = sample_name_match.group(1) if sample_name_match else "unknown_sample"

    # Define output chart filename
    output_filename = f"{sample_name}_blastn_top5hitspergene_species.png"

    # Plot pie chart
    fig, ax = plt.subplots()
    ax.pie(df['% of Total Hits'], labels=df['Species'], autopct='%1.1f%%', startangle=90)
    ax.set_title(f"Species Assignment – {sample_name}")
    plt.tight_layout()
    plt.savefig(output_filename, dpi=600)
    plt.close()

    log(f"✅ Saved pie chart as: {output_filename}")
    log(f"ℹ️ Hits: {len(species_hits)}, Empty: {empty_lines}, Malformed: {malformed_lines}, Unmatched: {unmatched_lines}")

log(f"\n🏁 All files processed. Log saved to: {log_file}\n")
