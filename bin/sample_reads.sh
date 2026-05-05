#!/usr/bin/env bash
# Randomly subsample FASTQ reads via staphb/seqtk:1.5 Docker container.
# For paired-end input, the same seed is applied to both files to preserve read pairing.
set -euo pipefail

NUM_READS=1000000
SEED=42
OUTDIR="data"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] R1.fastq.gz [R2.fastq.gz]

Randomly sample reads from FASTQ file(s) using seqtk (via Docker).
Paired-end: provide both R1 and R2; the same seed is used for both to preserve pairing.

Options:
  -n NUM   Number of reads to sample (default: ${NUM_READS})
  -s SEED  Random seed for reproducibility (default: ${SEED})
  -o DIR   Output directory (default: ${OUTDIR})
  -h       Show this help message
EOF
    exit 1
}

while getopts "n:s:o:h" opt; do
    case $opt in
        n) NUM_READS=$OPTARG ;;
        s) SEED=$OPTARG      ;;
        o) OUTDIR=$OPTARG    ;;
        h) usage             ;;
        *) usage             ;;
    esac
done
shift $((OPTIND - 1))

[[ $# -lt 1 ]] && usage

R1=$1
R2=${2:-}

command -v docker &>/dev/null || { echo "Error: docker not found"; exit 1; }

mkdir -p "$OUTDIR"
OUTDIR_ABS=$(realpath "$OUTDIR")

sample_file() {
    local infile=$1
    local outfile=$2
    local infile_abs indir inbase
    infile_abs=$(realpath "$infile")
    indir=$(dirname "$infile_abs")
    inbase=$(basename "$infile_abs")

    echo "  $(basename "$infile") -> $outfile"
    docker run --rm \
        --platform linux/amd64 \
        --user "$(id -u):$(id -g)" \
        -v "${indir}:/input:ro" \
        -v "${OUTDIR_ABS}:/output" \
        staphb/seqtk:1.5 \
        seqtk sample -s "$SEED" "/input/${inbase}" "$NUM_READS" \
        | gzip > "$outfile"
}

strip_ext() {
    local name
    name=$(basename "$1" .gz)
    name=$(basename "$name" .fastq)
    basename "$name" .fq
}

echo "Sampling ${NUM_READS} reads (seed=${SEED}) -> ${OUTDIR}/"
sample_file "$R1" "${OUTDIR_ABS}/$(strip_ext "$R1").fastq.gz"

if [[ -n "$R2" ]]; then
    sample_file "$R2" "${OUTDIR_ABS}/$(strip_ext "$R2").fastq.gz"
fi

echo "Done."
