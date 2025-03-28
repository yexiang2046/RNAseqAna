#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -g GEO_ACCESSION -o OUTPUT_DIR"
    echo "Download RNA-seq data from GEO/SRA"
    echo ""
    echo "Required arguments:"
    echo "  -g GEO_ACCESSION  GEO accession ID (e.g., GSE123456)"
    echo "  -o OUTPUT_DIR     Output directory for all files"
    echo "Optional arguments:"
    echo "  -h               Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hg:o:" opt; do
    case $opt in
        h) usage ;;
        g) GEO_ACC="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$GEO_ACC" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Check if required tools are installed
command -v fasterq-dump >/dev/null 2>&1 || { echo "Error: fasterq-dump is not installed. Please install SRA toolkit."; exit 1; }
command -v prefetch >/dev/null 2>&1 || { echo "Error: prefetch is not installed. Please install SRA toolkit."; exit 1; }

# Download GEO metadata
echo "Downloading GEO metadata..."
GEO_URL="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GEO_ACC}&targ=self&form=text&view=full"
curl -s "$GEO_URL" > "${GEO_ACC}_metadata.txt"

# Extract SRA IDs and sample information
echo "Extracting SRA IDs and sample information..."
awk '
    /^!Sample_geo_accession/ { split($0, a, "="); gsub(/^[ \t]+/, "", a[2]); sra[NR] = a[2] }
    /^!Sample_characteristics_ch1/ { split($0, a, "="); gsub(/^[ \t]+/, "", a[2]); info[NR] = a[2] }
    /^!Sample_title/ { split($0, a, "="); gsub(/^[ \t]+/, "", a[2]); title[NR] = a[2] }
    END {
        for (i in sra) {
            if (sra[i] != "") {
                print sra[i] "\t" title[i] "\t" info[i]
            }
        }
    }
' "${GEO_ACC}_metadata.txt" > sample_info.txt

# Create metadata.csv for edgeR
echo "Creating metadata.csv..."
echo "sample,condition" > metadata.csv
awk -F'\t' '
    {
        # Extract condition from characteristics (assuming format "condition: value")
        split($3, a, ":")
        if (length(a) > 1) {
            gsub(/^[ \t]+/, "", a[2])
            print $1 "," a[2]
        } else {
            # If no characteristics, use sample title
            print $1 "," $2
        }
    }
' sample_info.txt >> metadata.csv

# Create contrasts.csv for edgeR
echo "Creating contrasts.csv..."
echo "name,treatment,control" > contrasts.csv
# Get unique conditions
conditions=$(awk -F',' 'NR>1 {print $2}' metadata.csv | sort -u)
# Create contrasts between all pairs of conditions
first_condition=$(echo "$conditions" | head -n 1)
echo "$conditions" | tail -n +2 | while read condition; do
    echo "${condition}_vs_${first_condition},${condition},${first_condition}" >> contrasts.csv
done

# Download SRA data
echo "Downloading SRA data..."
mkdir -p raw_reads
cd raw_reads

# Download SRA data
while IFS=$'\t' read -r sra_id title info; do
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
done < ../sample_info.txt

cd ..

echo "Data download complete!"
echo "Files created:"
echo "1. ${GEO_ACC}_metadata.txt - Full GEO metadata"
echo "2. sample_info.txt - Extracted sample information"
echo "3. metadata.csv - Sample metadata for edgeR"
echo "4. contrasts.csv - Contrast definitions for edgeR"
echo "5. raw_reads/ - Raw FASTQ files"
echo ""
echo "You can now process the data using your preferred method." 