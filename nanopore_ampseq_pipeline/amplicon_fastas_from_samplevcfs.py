#!/usr/bin/env python3
"""
amplicon_fastas_from_samplevcfs.py

Generates per-sample and per-amplicon consensus FASTA files from a multi-sample VCF.

Steps per sample:
1. Extract the sample's variants from a multi-sample VCF
2. Filter to keep only variants where FMT/DP > 20
3. Normalize variants to split multiallelic and align indels
4. Generate a full consensus sequence, where missing data is output as N (using -M N)
5. Extract amplicon-specific sequences using BED regions
6. Write combined FASTA per amplicon across all samples

USAGE:
    python amplicon_fastas_from_samplevcfs.py \
        --bed regions.bed \
        --ref reference.fa \
        --vcf combined_filtered.vcf.gz \
        --outdir amplicons_merged_fastas

REQUIRES:
- bcftools
- samtools
- tabix

Sophie Moss, Campino Lab, Jun 2025
"""

import argparse
import subprocess
from pathlib import Path
from collections import defaultdict
from contextlib import redirect_stdout, redirect_stderr
import sys

class Tee:
    """Duplicate stdout/stderr to both terminal and a file."""
    def __init__(self, file, stream):
        self.file = file
        self.stream = stream

    def write(self, message):
        self.file.write(message)
        self.stream.write(message)

    def flush(self):
        self.file.flush()
        self.stream.flush()

def run(cmd, **kwargs):
    print(f"🔧 Running: {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, check=True, **kwargs)

def cleanup_files(*files):
    for f in files:
        p = Path(f)
        if p.exists():
            p.unlink()

def main():
    parser = argparse.ArgumentParser(description="Generate per-amplicon consensus FASTAs from a multi-sample VCF.")
    parser.add_argument("--bed", required=True, help="BED file with 4 columns: chrom, start, end, amplicon_name")
    parser.add_argument("--ref", required=True, help="Reference genome FASTA")
    parser.add_argument("--vcf", required=True, help="Combined multi-sample VCF (bgzipped)")
    parser.add_argument("--outdir", default="amplicons_merged_fastas", help="Directory for output FASTAs")
    args = parser.parse_args()

    bed = Path(args.bed)
    ref = Path(args.ref)
    vcf = Path(args.vcf)
    vcf_index = Path(str(vcf) + ".tbi")
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    # Redirect stdout and stderr to both terminal and log file
    log_file = outdir / "amplicons_merged_fastas.log"
    with open(log_file, "w") as log:
        tee = Tee(log, sys.stdout)
        with redirect_stdout(tee), redirect_stderr(tee):

            # Ensure reference is indexed
            if not Path(str(ref) + ".fai").exists():
                run(["samtools", "faidx", str(ref)])

            # Ensure VCF is indexed
            print("📦 Checking if VCF has index...")
            if not vcf_index.exists():
                print("🛠️  VCF index not found — creating VCF index")
                run(["tabix", "-p", "vcf", str(vcf)])
            else:
                print("✅ VCF index exists")

            # Get sample names from VCF
            result = subprocess.run(
                ["bcftools", "query", "-l", str(vcf)],
                capture_output=True, text=True, check=True
            )
            sample_names = result.stdout.strip().splitlines()

            # Read BED file
            amplicons = []
            with open(bed) as bed_file:
                for line in bed_file:
                    chrom, start, end, name = line.strip().split()
                    amplicons.append((chrom, start, end, name))

            # Track all sample sequences per amplicon
            combined_amplicon_seqs = defaultdict(list)

            for sample in sample_names:
                print(f"\n🧬 Processing sample: {sample}")
                sample_dir = outdir / sample
                sample_dir.mkdir(exist_ok=True)

                # STEP 1: Extract sample-specific VCF
                sample_vcf = sample_dir / f"{sample}.vcf.gz"
                run([
                    "bcftools", "view", "-s", sample, "-c", "1", str(vcf),
                    "-Oz", "-o", str(sample_vcf)
                ])
                run(["bcftools", "index", "-f", str(sample_vcf)])

                # STEP 2: Filter SNPs with FMT/DP > 20
                filtered_vcf = sample_dir / f"{sample}.FMTDP20.filtered.vcf.gz"
                run([
                    "bcftools", "view", "-i", "FMT/DP>20", "-v", "snps",
                    str(sample_vcf), "-Oz", "-o", str(filtered_vcf)
                ])
                run(["bcftools", "index", "-f", str(filtered_vcf)])

                # STEP 3: Normalize VCF
                norm_vcf = sample_dir / f"{sample}.FMTDP20.filtered.norm.vcf.gz"
                run([
                    "bcftools", "norm", "-m-any", "-f", str(ref),
                    str(filtered_vcf), "-Oz", "-o", str(norm_vcf)
                ])
                run(["bcftools", "index", "-f", str(norm_vcf)])

                # STEP 4: Generate consensus FASTA
                consensus_fa = sample_dir / f"{sample}_consensus.fa"
                print(f"🔧 Running: bcftools consensus -M N -f {ref} {norm_vcf}")
                with open(consensus_fa, "w") as out_fa:
                    subprocess.run([
                        "bcftools", "consensus", "-M", "N",
                        "-f", str(ref), str(norm_vcf)
                    ], check=True, stdout=out_fa)

                # STEP 5: Extract amplicon regions
                for chrom, start, end, name in amplicons:
                    region = f"{chrom}:{start}-{end}"
                    region_fa = sample_dir / f"{name}.fasta"

                    with open(region_fa, "w") as out_region:
                        run(["samtools", "faidx", str(consensus_fa), region], stdout=out_region)

                    seq_lines = Path(region_fa).read_text().splitlines()
                    seq = "".join([l for l in seq_lines if not l.startswith(">")])

                    if seq:
                        with open(region_fa, "w") as f:
                            f.write(f">{sample}\n{seq}\n")
                        combined_amplicon_seqs[name].append((sample, seq))
                        print(f"    ✔ {name} written for {sample}")
                    else:
                        print(f"    ⚠️  Empty sequence for {name} in {sample}, skipping.")

                cleanup_files(consensus_fa)

            # STEP 6: Write combined FASTAs per amplicon
            for name, entries in combined_amplicon_seqs.items():
                combined_path = outdir / f"{name}_combined.fasta"
                with open(combined_path, "w") as f:
                    for sample, seq in entries:
                        f.write(f">{sample}\n{seq}\n")
                print(f"📦 Combined FASTA written: {combined_path}")

            print(f"\n✅ All consensus sequences written to: {outdir}/")

if __name__ == "__main__":
    main()