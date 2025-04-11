#! /usr/bin/env python
import sys
import argparse
import subprocess as sp
import csv
import os
from uuid import uuid4

parser = argparse.ArgumentParser(description='Demultiplex plates')
parser.add_argument('-p', '--plate-layout', help='Input file', required=True)
parser.add_argument('-b', '--barcodes', help='Output file', required=True)
parser.add_argument('-f', '--fastq-dir', help='Fastq directory', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', default=4, type=int)
parser.add_argument('-o', '--output', default="samples.csv", help='Output file', required=True)
parser.add_argument('-m', '--max-mismatch', type=int, default=0, help='Maximum mismatches in barcode')

args = parser.parse_args()

# Read in plate layout
plate_layout = {}
for row in csv.DictReader(open(args.plate_layout, encoding='utf-8-sig')):
    plate_layout[row['Nanopore barcode ID']] = {k: v for k, v in row.items() if k not in ('Nanopore barcode ID', 'Well')}

# Read in barcodes
barcodes = {}
for row in csv.DictReader(open(args.barcodes, encoding='utf-8-sig')):
    barcodes[row['id']] = (row['forward'], row['reverse'])

# Demultiplex
output_rows = []
for bc in plate_layout:
    # Write temp barcode file
    tmp_barcode_file = f"{uuid4()}.csv"
    with open(tmp_barcode_file, 'w') as temp:
        rows = []
        for key, val in plate_layout[bc].items():
            new_id = f"{val}"
            rows.append({'id': new_id, 'forward': barcodes[key][0], 'reverse': barcodes[key][1]})
            output_rows.append({'sample': new_id, 'reads': new_id + ".fastq"})
        writer = csv.DictWriter(temp, fieldnames=['id', 'forward', 'reverse'])
        writer.writeheader()
        writer.writerows(rows)
    
    # ✅ Debugging print statement (INSIDE THE LOOP)
    print(f"\nTemporary barcode file {tmp_barcode_file} content:")
    with open(tmp_barcode_file, "r") as temp_file:
        print(temp_file.read())  # Print file contents for debugging

    # ✅ Check if Fastq file exists and is not empty
    fastq_file = f"{args.fastq_dir}/{bc}.fastq.gz"
    if not os.path.exists(fastq_file):
        print(f"❌ WARNING: Fastq file {fastq_file} is MISSING!")
        continue  # Skip processing this barcode
    elif os.path.getsize(fastq_file) == 0:
        print(f"⚠️ WARNING: Fastq file {fastq_file} is EMPTY!")
        continue  # Skip processing this barcode
    else:
        print(f"✅ Fastq file {fastq_file} exists and is ready for demultiplexing.")

    # Run demultiplexing
    sp.run(f"/mnt/storageG1/data/experiments/Exp149_eDNA_amplicon_Ghana/exp149/demux_nanopore_amplicon.py --max-mismatch {args.max_mismatch} --fastq {fastq_file} --barcodes {tmp_barcode_file} --edge-size 20 --log-prefix {bc}", shell=True)

    # Remove temp file
    sp.run(f"rm {tmp_barcode_file}", shell=True)

# Write final sample file
with open(args.output, "w") as O:
    writer = csv.DictWriter(O, fieldnames=['sample', 'reads'])
    writer.writeheader()
    writer.writerows(output_rows)