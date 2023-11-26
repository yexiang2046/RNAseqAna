#!/bin/bash

# script to count aligned reads in bam files using featureCounts

# Set Bash Options
set -Eeuo pipefail

# Function to display help message
display_help() {
    echo "Usage: $0 BAM_DIR OUT_FILE ANNOTATION_FILE OUT_DIR"
    echo "Arguments:"
    echo "  BAM_DIR    Path to the directory containing bam files"
    echo "  OUT_FILE    Output file name"
    echo "  ANNOTATION_FILE    Path to the annotation file"
    echo "  OUT_DIR    Path to the output directory"
    echo
}


# Function to log errors
log_error() {
    echo "Error: $1" >&2
}

# Trap errors and log them
trap 'log_error "Script error in line $LINENO"; exit 1' ERR

BAM_DIR=$1
OUT_FILE=$2
ANNOTATION_FILE=$3
OUT_DIR=$4  

# Check for the help option
if [[ "$1" == "--help" ]]; then
    display_help
    exit 0
fi


# count reads (paired-end)
featureCounts -T 16 -p -t exon -g gene_id -a ${ANNOTATION_FILE} -o ${OUT_DIR}/${OUT_FILE} ${BAM_DIR}/*.bam
