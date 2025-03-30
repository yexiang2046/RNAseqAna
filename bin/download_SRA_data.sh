#!/bin/bash

# Set strict error handling
set -e
set -o pipefail

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_file> -o <output_dir>"
    echo "  -i: Input file containing SRR IDs (one per line)"
    echo "  -o: Output directory for fastq files"
    exit 1
}

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Parse command line arguments
while getopts "i:o:" opt; do
    case $opt in
        i) input_file="$OPTARG";;
        o) output_dir="$OPTARG";;
        *) usage;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_file" ] || [ -z "$output_dir" ]; then
    usage
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    log "Error: Input file $input_file does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${output_dir}/raw_reads"

# Log start of processing
log "Starting SRA download process"
log "Input file: $input_file"
log "Output directory: $output_dir"

# Process each SRR ID
while read -r srr_id || [[ -n "$srr_id" ]]; do
    # Skip empty lines and comments
    [[ -z "$srr_id" || "$srr_id" =~ ^# ]] && continue
    
    log "Processing $srr_id"
    
    # Check if SRR ID is valid
    if [[ ! "$srr_id" =~ ^SRR[0-9]+$ ]]; then
        log "Warning: Invalid SRR ID format: $srr_id, skipping"
        continue
    fi
    
    # Download and convert to fastq
    log "Downloading $srr_id"
    if ! fasterq-dump "$srr_id" -O "${output_dir}/raw_reads" --split-files; then
        log "Error: Failed to download $srr_id"
        continue
    fi
    
    # Verify the download
    if ls "${output_dir}/raw_reads/${srr_id}"*.fastq 1> /dev/null 2>&1; then
        log "Successfully downloaded $srr_id"
    else
        log "Error: Output files not found for $srr_id"
    fi
done < "$input_file"

# Log completion
log "Download process completed"

# Final verification
fastq_count=$(ls "${output_dir}/raw_reads/"*.fastq 2>/dev/null | wc -l)
log "Total FASTQ files generated: $fastq_count"

if [ "$fastq_count" -eq 0 ]; then
    log "Error: No FASTQ files were generated"
    exit 1
fi

exit 0 