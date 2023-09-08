#!/bin/bash

# Set Bash options for error handling and command execution 
set -Eeuo pipefail

# Function to display help message
display_help() {
    echo "Usage: $0 FILE_PATH READ1_PAT READ2_PAT"
    echo "Arguments:"
    echo "  FILE_PATH    Path to the directory containing fastq.gz files"
    echo "  READ1_PAT    Pattern for Read 1 fastq.gz files"
    echo "  READ2_PAT    Pattern for Read 2 fastq.gz files"
    echo
}

# Function to log errors
log_error() {
    echo "Error: $1" >&2
}

# Trap errors and log them
trap 'log_error "Script error in line $LINENO"; exit 1' ERR


# Check for the help option
if [[ "$1" == "--help" ]]; then
    display_help
    exit 0
fi



# assign input variables
FILE_PATH=$1

READ1_PAT=$2
READ2_PAT=$3



samples_id=$(find "$FILE_PATH" -maxdepth 1 -type f -name '*.fastq.gz' -printf '%P\n' | awk -F'_' '{print $1}' | sort | uniq)
printf '%s\n' "$samples_id" | while IFS= read -r f; do
    echo "Read1: ${FILE_PATH}/$f${READ1_PAT}"
    echo "Read2: ${FILE_PATH}/$f${READ2_PAT}"
    fastp -w 16 -l 20 -i ${FILE_PATH}/$f${READ1_PAT} -I ${FILE_PATH}/$f${READ2_PAT} -o ${FILE_PATH}/$f.R1.fastp.fastq.gz -O ${FILE_PATH}/$f.R2.fastp.fastq.gz
done


