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

# Download SRA metadata from Run Selector
SRA_URL="https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${GEO_ACC}&output=file"
curl -s "$SRA_URL" > "${GEO_ACC}_metadata.txt"

# Extract SRA IDs and sample information
awk -F',' '
    NR==1 {
        for(i=1; i<=NF; i++) {
            if($i=="Run") run_col=i
            if($i=="Sample Name") sample_col=i
            if($i=="condition") condition_col=i
            if($i=="cell_type") cell_type_col=i
            if($i=="tissue") tissue_col=i
            if($i=="time_point") time_point_col=i
        }
    }
    NR>1 {
        # Clean up the fields (remove quotes and escape characters)
        gsub(/"/, "", $run_col)
        gsub(/"/, "", $sample_col)
        gsub(/"/, "", $condition_col)
        gsub(/"/, "", $cell_type_col)
        gsub(/"/, "", $tissue_col)
        gsub(/"/, "", $time_point_col)
        gsub(/,/, "_", $sample_col)
        gsub(/,/, "_", $condition_col)
        gsub(/,/, "_", $cell_type_col)
        gsub(/,/, "_", $tissue_col)
        gsub(/,/, "_", $time_point_col)
        
        # Combine condition, cell_type, tissue, and time_point for a more descriptive condition
        condition = $condition_col
        if ($cell_type_col != "") condition = condition "_" $cell_type_col
        if ($tissue_col != "") condition = condition "_" $tissue_col
        if ($time_point_col != "") condition = condition "_" $time_point_col
        
        print $run_col "\t" $sample_col "\t" condition
    }
' "${GEO_ACC}_metadata.txt" > sample_info.txt

# Create metadata.csv for edgeR
echo "sample,condition" > metadata.csv
awk -F'\t' '
    {
        print $1 "," $3
    }
' sample_info.txt >> metadata.csv

# Create contrasts.csv for edgeR
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
echo "1. ${GEO_ACC}_metadata.txt - Full SRA metadata"
echo "2. sample_info.txt - Extracted sample information"
echo "3. metadata.csv - Sample metadata for edgeR"
echo "4. contrasts.csv - Contrast definitions for edgeR"
echo "5. raw_reads/ - Raw FASTQ files"
echo ""
echo "You can now process the data using your preferred method." 