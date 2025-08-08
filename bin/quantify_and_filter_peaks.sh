#!/usr/bin/env bash

# Peak Quantification and Filtering Script
# This script quantifies peaks from IDR output and filters them based on read counts

set -euo pipefail

# Help message
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -i IDR_PEAKS -b BAM_FILE -o OUTPUT_DIR

Quantify peaks from IDR output and filter based on read counts.

Required arguments:
  -i IDR_PEAKS    Input IDR peaks file (BED format)
  -b BAM_FILE     Input BAM file for read counting
  -o OUTPUT_DIR   Output directory for results

Optional arguments:
  -m MIN_READS    Minimum read count per peak (default: 10)
  -w MIN_WIDTH    Minimum peak width in bp (default: 50)
  -W MAX_WIDTH    Maximum peak width in bp (default: 10000)
  -g GTF_FILE     GTF file for gene annotation (optional)
  -s SAMPLE_ID    Sample ID for output files (default: extracted from BAM filename)
  -h              Show this help message

Output files:
  - quantified_peaks.bed: Peaks with read counts and statistics
  - filtered_peaks.bed: Peaks passing filtering criteria
  - peak_stats.txt: Summary statistics
  - read_counts_per_peak.txt: Detailed read count information
  - annotated_peaks.bed: Peaks with gene annotations (if GTF provided)
  - peak_gene_association.txt: Peak-gene associations (if GTF provided)

Example:
  $0 -i idr_peaks.bed -b sample.bam -o results/ -m 20 -w 100 -g genes.gtf
EOF
    exit 1
}

# Parse command line arguments
while getopts "hi:b:o:m:w:W:g:s:" opt; do
    case $opt in
        h) usage ;;
        i) IDR_PEAKS="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) MIN_READS="$OPTARG" ;;
        w) MIN_WIDTH="$OPTARG" ;;
        W) MAX_WIDTH="$OPTARG" ;;
        g) GTF_FILE="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "${IDR_PEAKS:-}" ] || [ -z "${BAM_FILE:-}" ] || [ -z "${OUTPUT_DIR:-}" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input files exist
if [ ! -f "$IDR_PEAKS" ]; then
    echo "Error: IDR peaks file does not exist: $IDR_PEAKS"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file does not exist: $BAM_FILE"
    exit 1
fi

# Set default values
MIN_READS=${MIN_READS:-10}
MIN_WIDTH=${MIN_WIDTH:-50}
MAX_WIDTH=${MAX_WIDTH:-10000}
SAMPLE_ID=${SAMPLE_ID:-$(basename "$BAM_FILE" .bam)}

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create temporary directory
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

echo "=== Peak Quantification and Filtering ==="
echo "Sample ID: $SAMPLE_ID"
echo "Input IDR peaks: $IDR_PEAKS"
echo "Input BAM file: $BAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Minimum reads per peak: $MIN_READS"
echo "Peak width range: ${MIN_WIDTH}-${MAX_WIDTH} bp"
echo ""

# Step 1: Convert BAM to BED for read counting
echo "Step 1: Converting BAM to BED format..."
bedtools bamtobed -i "$BAM_FILE" > "$TEMP_DIR/reads.bed"

# Step 2: Count reads overlapping with each peak
echo "Step 2: Counting reads per peak..."
bedtools coverage -a "$IDR_PEAKS" -b "$TEMP_DIR/reads.bed" -counts > "$TEMP_DIR/peaks_with_counts.bed"

# Step 3: Add peak information and calculate statistics
echo "Step 3: Calculating peak statistics..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    # Calculate peak width
    width = $3 - $2
    
    # Calculate peak center
    center = int(($2 + $3) / 2)
    
    # Calculate peak density (reads per bp)
    density = $NF / width
    
    # Print: chr, start, end, peak_name, score, strand, read_count, width, center, density
    print $1, $2, $3, $4, $5, $6, $NF, width, center, density
}' "$TEMP_DIR/peaks_with_counts.bed" > "$OUTPUT_DIR/${SAMPLE_ID}_quantified_peaks.bed"

# Step 4: Filter peaks based on criteria
echo "Step 4: Filtering peaks based on read count and width criteria..."
awk -F'\t' -v min_reads="$MIN_READS" -v min_width="$MIN_WIDTH" -v max_width="$MAX_WIDTH" '
BEGIN {OFS="\t"}
{
    read_count = $7
    width = $8
    
    # Filter criteria: minimum reads, width constraints
    if (read_count >= min_reads && width >= min_width && width <= max_width) {
        print $0
    }
}' "$OUTPUT_DIR/${SAMPLE_ID}_quantified_peaks.bed" > "$OUTPUT_DIR/${SAMPLE_ID}_filtered_peaks.bed"

# Step 5: Generate statistics
echo "Step 5: Generating peak statistics..."
{
    echo "=== Peak Quantification Statistics for $SAMPLE_ID ==="
    echo "Input peaks: $(wc -l < "$IDR_PEAKS")"
    echo "Peaks with read counts: $(wc -l < "$OUTPUT_DIR/${SAMPLE_ID}_quantified_peaks.bed")"
    echo "Filtered peaks: $(wc -l < "$OUTPUT_DIR/${SAMPLE_ID}_filtered_peaks.bed")"
    echo ""
    echo "=== Read Count Statistics ==="
    awk -F'\t' '{print $7}' "$OUTPUT_DIR/${SAMPLE_ID}_quantified_peaks.bed" | sort -n | awk '
    {
        count[++n] = $1
        sum += $1
    }
    END {
        if (n > 0) {
            mean = sum / n
            if (n % 2 == 1) {
                median = count[(n + 1) / 2]
            } else {
                median = (count[n/2] + count[n/2 + 1]) / 2
            }
            print "Mean read count per peak: " mean
            print "Median read count per peak: " median
            print "Min read count: " count[1]
            print "Max read count: " count[n]
            print "Total reads in peaks: " sum
        }
    }'
    echo ""
    echo "=== Peak Width Statistics ==="
    awk -F'\t' '{print $8}' "$OUTPUT_DIR/${SAMPLE_ID}_quantified_peaks.bed" | sort -n | awk '
    {
        count[++n] = $1
        sum += $1
    }
    END {
        if (n > 0) {
            mean = sum / n
            if (n % 2 == 1) {
                median = count[(n + 1) / 2]
            } else {
                median = (count[n/2] + count[n/2 + 1]) / 2
            }
            print "Mean peak width: " mean " bp"
            print "Median peak width: " median " bp"
            print "Min peak width: " count[1] " bp"
            print "Max peak width: " count[n] " bp"
        }
    }'
    echo ""
    echo "=== Filtering Criteria ==="
    echo "Minimum reads per peak: $MIN_READS"
    echo "Peak width range: ${MIN_WIDTH}-${MAX_WIDTH} bp"
} > "$OUTPUT_DIR/${SAMPLE_ID}_peak_stats.txt"

# Step 6: Create detailed read count file
echo "Step 6: Creating detailed read count file..."
{
    echo -e "peak_id\tchrom\tstart\tend\tread_count\tpeak_width\tpeak_center\tdensity"
    awk -F'\t' 'BEGIN {OFS="\t"} {
        print $4, $1, $2, $3, $7, $8, $9, $10
    }' "$OUTPUT_DIR/${SAMPLE_ID}_quantified_peaks.bed"
} > "$OUTPUT_DIR/${SAMPLE_ID}_read_counts_per_peak.txt"

# Step 7: Gene annotation (if GTF file provided)
if [ -n "${GTF_FILE:-}" ]; then
    if [ -f "$GTF_FILE" ]; then
        echo "Step 7: Annotating peaks with genes..."
        
        # Convert GTF to BED format for annotation
        awk -F'\t' '$3 == "gene" {
            # Extract gene ID from attributes
            split($9, attrs, ";")
            gene_id = ""
            gene_name = ""
            for (i in attrs) {
                if (attrs[i] ~ /gene_id/) {
                    gsub(/gene_id "/, "", attrs[i])
                    gsub(/"/, "", attrs[i])
                    gene_id = attrs[i]
                }
                if (attrs[i] ~ /gene_name/) {
                    gsub(/gene_name "/, "", attrs[i])
                    gsub(/"/, "", attrs[i])
                    gene_name = attrs[i]
                }
            }
            print $1 "\t" $4-1 "\t" $5 "\t" gene_id "\t" gene_name
        }' "$GTF_FILE" > "$TEMP_DIR/genes.bed"
        
        # Find closest genes to peaks
        bedtools closest -a "$OUTPUT_DIR/${SAMPLE_ID}_filtered_peaks.bed" -b "$TEMP_DIR/genes.bed" -d > "$TEMP_DIR/peaks_with_genes.bed"
        
        # Create annotated peaks file
        awk -F'\t' 'BEGIN {OFS="\t"} {
            # Format: chr, start, end, peak_id, score, strand, read_count, width, center, density, gene_id, gene_name, distance
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $15, $16
        }' "$TEMP_DIR/peaks_with_genes.bed" > "$OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks.bed"
        
        # Create peak-gene association file
        {
            echo -e "peak_id\tchrom\tpeak_start\tpeak_end\tread_count\tgene_id\tgene_name\tdistance_to_gene"
            awk -F'\t' 'BEGIN {OFS="\t"} {
                print $4, $1, $2, $3, $7, $11, $12, $13
            }' "$TEMP_DIR/peaks_with_genes.bed"
        } > "$OUTPUT_DIR/${SAMPLE_ID}_peak_gene_association.txt"
        
        echo "Gene annotation completed."
    else
        echo "Warning: GTF file does not exist: $GTF_FILE"
    fi
fi

echo ""
echo "=== Summary ==="
echo "Processing completed for sample: $SAMPLE_ID"
echo "Output files in: $OUTPUT_DIR"
echo ""
echo "Generated files:"
echo "  - ${SAMPLE_ID}_quantified_peaks.bed: All peaks with read counts"
echo "  - ${SAMPLE_ID}_filtered_peaks.bed: Filtered peaks"
echo "  - ${SAMPLE_ID}_peak_stats.txt: Summary statistics"
echo "  - ${SAMPLE_ID}_read_counts_per_peak.txt: Detailed read counts"
if [ -n "${GTF_FILE:-}" ] && [ -f "$GTF_FILE" ]; then
    echo "  - ${SAMPLE_ID}_annotated_peaks.bed: Peaks with gene annotations"
    echo "  - ${SAMPLE_ID}_peak_gene_association.txt: Peak-gene associations"
fi
echo ""
echo "Filtering criteria applied:"
echo "  - Minimum reads per peak: $MIN_READS"
echo "  - Peak width range: ${MIN_WIDTH}-${MAX_WIDTH} bp"
echo ""
echo "Peak quantification and filtering completed successfully!" 