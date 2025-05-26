#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i INPUT_TSV -o OUTPUT_BED [-s SAMPLE_NAME]"
    echo "Transform Piranha peak calling output (TSV) to BED format"
    echo ""
    echo "Arguments:"
    echo "  -i INPUT_TSV    Input TSV file from Piranha"
    echo "  -o OUTPUT_BED   Output BED file"
    echo "  -s SAMPLE_NAME  Sample name (default: derived from input filename)"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:o:s:" opt; do
    case $opt in
        h) usage ;;
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input file exists
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file does not exist: $INPUT"
    exit 1
fi

# If sample name not provided, derive from input filename
if [ -z "$SAMPLE" ]; then
    SAMPLE=$(basename "$INPUT" .tsv | sed 's/_peaks$//')
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT")
mkdir -p "$OUTPUT_DIR"

echo "Converting Piranha peaks to BED format..."
echo "Input: $INPUT"
echo "Output: $OUTPUT"
echo "Sample: $SAMPLE"

# Convert TSV to BED format
# Expected Piranha TSV format: chr start end read_count total_count strand p-value
# Output BED format: chr start end name score strand
awk -v OFS='\t' -v sample="$SAMPLE" '
    # Skip header
    NR>1 {
        chr=$1
        start=$2
        end=$3
        read_count=$4
        strand=$6
        pval=$7
        
        # Format name as sample_peak_number_pvalue
        name=sprintf("%s_peak_%d_pval%.2e", sample, NR-1, pval)
        
        # Use read_count as score
        score=read_count
        
        # Print BED format
        print chr, start, end, name, score, strand
    }' "$INPUT" > "$OUTPUT"

# Count peaks
PEAK_COUNT=$(wc -l < "$OUTPUT")
echo "Converted $PEAK_COUNT peaks"

# Basic validation
echo "Validating output..."
if [ "$PEAK_COUNT" -eq 0 ]; then
    echo "Warning: No peaks found in output file"
    exit 1
fi

# Print first few peaks as example
echo "First few peaks:"
head -n 3 "$OUTPUT"

echo "Done! Output written to: $OUTPUT" 