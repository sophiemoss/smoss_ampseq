#!/bin/bash
set -e

# === Configuration ===
DB_PATH="/mnt/storage13/vectors/DenovoAssemblyTools/runs/realone/nt/nt"
THREADS_PER_JOB=20
FASTA_DIR="."
N_CORES=$(nproc)
MAX_PROCS=$((N_CORES / THREADS_PER_JOB))

echo "🔍 Starting BLAST using xargs with $MAX_PROCS parallel jobs..."

run_blast() {
    file="$1"
    base=$(basename "$file" .fasta)

    echo "🚀 Running BLAST for $base"

    blastn -query "$file" \
        -db "$DB_PATH" \
        -evalue 1e-10 \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -num_threads "$THREADS_PER_JOB" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -out "${base}_nt_results.txt"

    echo "✅ Finished $base"
}

export -f run_blast
export DB_PATH THREADS_PER_JOB

find "$FASTA_DIR" -name "*.fasta" | xargs -n 1 -P "$MAX_PROCS" -I {} bash -c 'run_blast "$@"' _ {}

echo "🎯 All BLAST jobs completed (xargs)."
