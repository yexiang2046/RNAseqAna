#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i ID_LIST -o OUTPUT_DIR"
    echo "Download RNA-seq data from SRA"
    echo ""
    echo "Required arguments:"
    echo "  -i ID_LIST     Input file containing SRA IDs (one per line)"
    echo "  -o OUTPUT_DIR  Output directory for all files"
    echo "Optional arguments:"
    echo "  -h            Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:o:" opt; do
    case $opt in
        h) usage ;;
        i) ID_LIST="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$ID_LIST" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

if [ ! -f "$ID_LIST" ]; then
    echo "Error: ID_LIST file does not exist"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Check if required tools are installed
command -v fasterq-dump >/dev/null 2>&1 || { echo "Error: fasterq-dump is not installed. Please install SRA toolkit."; exit 1; }
command -v prefetch >/dev/null 2>&1 || { echo "Error: prefetch is not installed. Please install SRA toolkit."; exit 1; }

# Download SRA data
echo "Downloading SRA data..."
mkdir -p "${OUTPUT_DIR}/raw_reads"
cd "${OUTPUT_DIR}/raw_reads"

# Download SRA data
while read -r sra_id; do
    if [ ! -z "$sra_id" ]; then
        echo "Downloading $sra_id..."
        # First prefetch the SRA data
        prefetch "$sra_id"
        
        # Then use fasterq-dump with common parameters
        fasterq-dump \
            --split-files \
            --threads 4 \
            --outdir . \
            --skip-technical \
            --read-filter pass \
            --min-read-len 50 \
            --max-read-len 1000 \
            --include-technical \
            --force \
            "$sra_id"
    fi
done < ../../"$ID_LIST"

cd ..

echo "Data download complete!"
echo "Files created:"
echo "1. raw_reads/ - Raw FASTQ files"
echo ""
echo "You can now process the data using your preferred method." 