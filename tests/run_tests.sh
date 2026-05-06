#!/usr/bin/env bash
set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

STAR_INDEX=""
GTF=""
DATA_DIR="tests_data"
OUTDIR="tests_results"

usage() {
    cat <<EOF
Usage: $(basename "$0") --star_index DIR --gtf FILE [--data_dir DIR] [--outdir DIR]

Run the RNA-seq pipeline against sampled test FASTQs and validate outputs.

Required:
  --star_index DIR   Path to a pre-built STAR index directory
  --gtf FILE         Path to the annotation GTF file

Optional:
  --data_dir  DIR    Directory containing paired-end FASTQ files (default: tests_data)
  --outdir    DIR    Pipeline output directory (default: tests_results)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --star_index) STAR_INDEX=$2; shift 2 ;;
        --gtf)        GTF=$2;        shift 2 ;;
        --data_dir)   DATA_DIR=$2;   shift 2 ;;
        --outdir)     OUTDIR=$2;     shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

[[ -z "$STAR_INDEX" || -z "$GTF" ]] && usage

PASS=0
FAIL=0

check() {
    local label=$1
    local pattern=$2
    if compgen -G "$pattern" > /dev/null 2>&1; then
        echo -e "  ${GREEN}PASS${NC}  $label"
        ((PASS++))
    else
        echo -e "  ${RED}FAIL${NC}  $label  ($pattern not found)"
        ((FAIL++))
    fi
}

echo "Running pipeline..."
echo "  data_dir:   $DATA_DIR"
echo "  star_index: $STAR_INDEX"
echo "  gtf:        $GTF"
echo "  outdir:     $OUTDIR"
echo ""

nextflow run "$(dirname "$0")/../main.nf" \
    --data_dir   "$DATA_DIR"   \
    --star_index "$STAR_INDEX" \
    --gtf        "$GTF"        \
    --outdir     "$OUTDIR"     \
    -resume

echo ""
echo "Validating outputs..."
check "trimmed reads"     "${OUTDIR}/trimmed/*.fastq.gz"
check "aligned BAM"       "${OUTDIR}/aligned/*.bam"
check "feature counts"    "${OUTDIR}/feature_counts/counts.txt"
check "MultiQC report"    "${OUTDIR}/multiqc_report/multiqc_report.html"

echo ""
if [[ $FAIL -eq 0 ]]; then
    echo -e "${GREEN}All ${PASS} checks passed.${NC}"
    exit 0
else
    echo -e "${RED}${FAIL} check(s) failed.${NC}"
    exit 1
fi
