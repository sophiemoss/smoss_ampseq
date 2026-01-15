#!/bin/bash
set -euo pipefail

# Pipeline: ONT variant calling and consensus sequence generation using Clair3
#
# Outputs layout:
#   sample_consensus_amplicon_fastas/
#     ├── <sample>/
#     │     ├── <sample>.<region>.hap1.fa
#     │     └── <sample>.<region>.hap2.fa
#     └── combined.<region>.hap1.fasta
#     └── combined.<region>.hap2.fasta
#
# Required flags:
#   -r, --reference  <FASTA>
#   -g, --regions    <BED>
#
# Optional env vars:
#   MASK_DP (default 5), THREADS (default 8), CLAIR3_MODEL (see default below)
#
# Example run: bash /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/consensus_fastas_loop.sh -r Anopheles_gambiae.AgamP4.dna.toplevel.fa -g regions.bed 2>&1 > log.txt

# ---------- CLI parsing ----------
REFERENCE=""
REGIONS_BED=""
THREADS=${THREADS:-8}
CLAIR3_MODEL=${CLAIR3_MODEL:-/mnt/storage11/sophie/miniconda3/envs/ont/bin/models/r1041_e82_400bps_sup_v500}
MASK_DP=${MASK_DP:-5}
OUTDIR="sample_consensus_amplicon_fastas"

usage() {
  echo "Usage: $0 -r <reference.fa> -g <regions.bed>"
  echo "  -r, --reference   Path to reference FASTA"
  echo "  -g, --regions     Path to regions BED (BED3/BED4)"
  echo "  [env] MASK_DP     Low-depth threshold for masking (default: ${MASK_DP})"
  echo "  [env] THREADS     Threads for Clair3 (default: ${THREADS})"
  echo "  [env] CLAIR3_MODEL  Clair3 model path (default: ${CLAIR3_MODEL})"
  exit 1
}

# Parse long/short flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reference) REFERENCE="$2"; shift 2 ;;
    -g|--regions)   REGIONS_BED="$2"; shift 2 ;;
    -h|--help)      usage ;;
    --) shift; break ;;
    -*) echo "ERROR: Unknown option: $1"; usage ;;
    *)  echo "ERROR: Unexpected positional argument: $1"; usage ;;
  esac
done

# ---------- Main loop  ----------
for bam in *.bam; do
  [ -e "$bam" ] || { echo "No BAM files found in $(pwd)"; break; }
  sample=$(basename "$bam" .bam)
  sample_dir="${OUTDIR}/${sample}"
  mkdir -p "$sample_dir"

  (
    echo "🧬 Processing sample: $sample"

    # Ensure BAM index exists
    [ -f "${bam}.bai" ] || samtools index "$bam"

    # 1) Variant calling (Clair3)
    run_clair3.sh \
      --bam_fn "$bam" \
      --ref_fn "$REFERENCE" \
      --threads "$THREADS" \
      --platform ont \
      --output "${sample_dir}/${sample}.clair3" \
      --model_path "$CLAIR3_MODEL" \
      --include_all_ctgs \
      --bed_fn="$REGIONS_BED" \
      --sample_name="$sample"

    # 2) Rename Clair3 outputs to include sample name (skip if missing)
    for f in merge_output full_alignment pileup; do
      src_vcf="${sample_dir}/${sample}.clair3/${f}.vcf.gz"
      src_tbi="${src_vcf}.tbi"
      [ -f "$src_vcf" ] && mv "$src_vcf" "${sample_dir}/${sample}.${f}.vcf.gz"
      [ -f "$src_tbi" ] && mv "$src_tbi" "${sample_dir}/${sample}.${f}.vcf.gz.tbi"
    done

    # 3) Pick VCF: always use Clair3's merge_output; skip sample if missing/empty
    in_vcf="${sample_dir}/${sample}.merge_output.vcf.gz"
    if [ ! -s "$in_vcf" ]; then
      echo "⚠️  $sample: Clair3 merge_output VCF missing or empty (${in_vcf}). Skipping."
      exit 0
    fi
    tabix -f "$in_vcf" 2>/dev/null || bcftools index -t -f "$in_vcf"

    # Determine sample name from VCF (single-sample expected)
    vcf_sample=$(bcftools query -l "$in_vcf" | head -1 || true)

    # 4) Phase variants (NO normalization; phase directly on merge_output VCF)
    #    If Whatshap fails, fall back to unphased VCF.
    phased_vcf="${sample_dir}/${sample}.phased.vcf.gz"
    if [ -n "$vcf_sample" ]; then
      whatshap phase \
        --reference "$REFERENCE" \
        --indels \
        --sample "$vcf_sample" \
        --ignore-read-groups \
        -o "$phased_vcf" \
        "$in_vcf" "$bam" \
      || { echo "⚠️  $sample: Whatshap failed; using unphased VCF for consensus."; cp -f "$in_vcf" "$phased_vcf"; }
    else
      echo "⚠️  $sample: Could not determine sample name from VCF; skipping phasing."
      cp -f "$in_vcf" "$phased_vcf"
    fi
    tabix -f "$phased_vcf" 2>/dev/null || bcftools index -t -f "$phased_vcf"

    # 5) BAM-based low-depth mask
    {
      T="$MASK_DP"
      samtools depth -a "$bam" \
      | awk -v T="$T" 'BEGIN{OFS="\t"}
        {
          if ($3<T) {
            if (s=="") { chr=$1; s=$2-1; e=$2 }
            else if ($1==chr && $2==e+1) { e=$2 }
            else { print chr,s,e; chr=$1; s=$2-1; e=$2 }
          } else if (s!="") { print chr,s,e; s="" }
        }
        END{ if (s!="") print chr,s,e }'
    } > "${sample_dir}/${sample}.depthmask.bed"

    # 6) VCF no-call mask to force any positions with ./. or .|. to be N in the consensus

    if [ -n "${vcf_sample:-}" ]; then
      bcftools query -s "$vcf_sample" \
        -i 'GT=".|." || GT="./."' \
        -f '%CHROM\t%POS0\t%POS\n' "$phased_vcf" \
        > "${sample_dir}/${sample}.nocall.mask.bed" || true
    else
      : > "${sample_dir}/${sample}.nocall.mask.bed"
    fi

    # 7) Merge depth + no-call masks to create one mask with positions which will be forced to N
    if command -v bedtools >/dev/null 2>&1; then
      cat "${sample_dir}/${sample}.depthmask.bed" "${sample_dir}/${sample}.nocall.mask.bed" \
        | sort -k1,1 -k2,2n \
        | bedtools merge > "${sample_dir}/${sample}.mask.bed"
    else
      cat "${sample_dir}/${sample}.depthmask.bed" "${sample_dir}/${sample}.nocall.mask.bed" \
        | sort -k1,1 -k2,2n \
        | awk 'BEGIN{OFS="\t"}
            NR==1 {c=$1; s=$2; e=$3; next}
            { if($1==c && $2<=e){ if($3>e) e=$3 }
              else { print c,s,e; c=$1; s=$2; e=$3 } }
            END{ if(NR>0) print c,s,e }' > "${sample_dir}/${sample}.mask.bed"
    fi

   # 8) Build per-amplicon haplotype consensus files
    # ensure sample_arg is defined earlier:
    # vcf_sample=$(bcftools query -l "$in_vcf" | head -1 || true)
    # sample_arg=(); [ -n "${vcf_sample:-}" ] && sample_arg=(-s "$vcf_sample")

    while IFS=$'\t' read -r chrom start end name rest; do
      region="${chrom}:${start}-${end}"

      # skip if region not in reference
      samtools faidx "$REFERENCE" "$region" >/dev/null 2>&1 || {
        echo "⚠️  $sample: region $region not in reference; skipping."
        continue
      }

      for H in 1 2; do
        out="${sample_dir}/${sample}.${name}.hap${H}.fa"
        samtools faidx "$REFERENCE" "$region" \
        | bcftools consensus \
            "${sample_arg[@]}" \
            -H "$H" \
            -m "${sample_dir}/${sample}.mask.bed" \
            "$phased_vcf" \
        | tr '*' 'N' \
        | sed "1s/^>.*/>${sample}_${name}_hap${H}/" \
        > "$out"

        [ -s "$out" ] || echo "⚠️  $sample: empty consensus for $region (hap${H})."
      done
    done < "$REGIONS_BED"

    echo "✅ Finished processing $sample"
  ) || { echo "⚠️  Sample '$sample' failed — continuing with the next sample."; continue; }

done

# 10) Concatenate per-amplicon FASTAs (hap1/hap2)
echo "📦 Concatenating per-amplicon FASTAs..."

# For files: $OUTDIR/<sample>/<sample>.<amplicon>.<hap>.fa
find "$OUTDIR" -type f -name '*.fa' | while read -r f; do
  sample=$(basename "$(dirname "$f")")   # sample from parent dir
  base=$(basename "$f" .fa)              # sample.amplicon.hap
  rest=${base#"$sample."}                # amplicon.hap
  hap=${rest##*.}                        # hap1 | hap2
  amp=${rest%.*}                         # amplicon (may contain dots)
  out="${OUTDIR}/combined.${amp}.${hap}.fasta"
  sed "1s/^>.*/>${sample}_${amp}_${hap}/" "$f" >> "$out"
done

echo "✅ Combined FASTAs written to $OUTDIR (combined.<amp>.<hap>.fasta)."

echo "✅ Combined FASTA files created in ${OUTDIR}:"
find "$OUTDIR" -maxdepth 1 -type f -name "combined.*.fasta" -print
