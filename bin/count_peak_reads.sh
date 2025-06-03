#!/bin/bash

# Usage: count_peak_reads.sh -b BAM_FILE -p PEAKS_FILE -c CHROMOSOME -o OUTPUT_PREFIX

# Parse command line arguments
while getopts "b:p:c:o:" opt; do
    case $opt in
        b) BAM_FILE="$OPTARG";;
        p) PEAKS_FILE="$OPTARG";;
        c) CHROMOSOME="$OPTARG";;
        o) OUTPUT_PREFIX="$OPTARG";;
        \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
    esac
done

# Check required arguments
if [ -z "$BAM_FILE" ] || [ -z "$PEAKS_FILE" ] || [ -z "$CHROMOSOME" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Usage: $0 -b BAM_FILE -p PEAKS_FILE -c CHROMOSOME -o OUTPUT_PREFIX"
    exit 1
fi

# Function to check if a file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 does not exist"
        exit 1
    fi
}

# Check input files
check_file "$BAM_FILE"
check_file "$PEAKS_FILE"

# Sort and index BAM file if not already indexed
if [ ! -f "${BAM_FILE}.bai" ]; then
    echo "Indexing BAM file..."
    samtools sort -@ 4 -o "${OUTPUT_PREFIX}_sorted.bam" "$BAM_FILE"
    samtools index -@ 4 "${OUTPUT_PREFIX}_sorted.bam"
    BAM_FILE="${OUTPUT_PREFIX}_sorted.bam"
fi

# Sort peaks file if not already sorted
if [ ! -f "${PEAKS_FILE}.sorted" ]; then
    echo "Sorting peaks file..."
    sort -k1,1 -k2,2n "$PEAKS_FILE" > "${PEAKS_FILE}.sorted"
    PEAKS_FILE="${PEAKS_FILE}.sorted"
fi

# Count reads in peaks using bedtools
echo "Counting reads in peaks..."
bedtools coverage -a "$PEAKS_FILE" -b "$BAM_FILE" -counts -S > "${OUTPUT_PREFIX}_peak_counts.txt"

# Calculate viral genome coverage
echo "Calculating viral genome coverage..."
samtools depth -r "$CHROMOSOME" "$BAM_FILE" > "${OUTPUT_PREFIX}_viral_depth.txt"

# Get viral genome length
VIRAL_LENGTH=$(samtools view -H "$BAM_FILE" | grep "SN:${CHROMOSOME}" | awk '{print $3}' | sed 's/LN://')

# Calculate coverage statistics
echo "Calculating coverage statistics..."
awk -v chr="$CHROMOSOME" -v len="$VIRAL_LENGTH" '
BEGIN {
    total_bases = 0
    covered_bases = 0
    total_depth = 0
    max_depth = 0
}
{
    if ($3 > 0) {
        covered_bases++
        total_depth += $3
        if ($3 > max_depth) max_depth = $3
    }
    total_bases++
}
END {
    avg_depth = (total_bases > 0) ? total_depth/total_bases : 0
    coverage = (total_bases > 0) ? (covered_bases/total_bases)*100 : 0
    print "Viral Genome Coverage Statistics"
    print "================================"
    print "Chromosome: " chr
    print "Genome length: " len " bp"
    print "Total bases sequenced: " total_bases
    print "Bases with coverage: " covered_bases
    print "Average depth: " avg_depth
    print "Maximum depth: " max_depth
    print "Coverage percentage: " coverage "%"
    print "Total reads: " total_depth
}' "${OUTPUT_PREFIX}_viral_depth.txt" > "${OUTPUT_PREFIX}_viral_coverage.txt"

# Clean up temporary files
if [ -f "${OUTPUT_PREFIX}_sorted.bam" ]; then
    rm "${OUTPUT_PREFIX}_sorted.bam" "${OUTPUT_PREFIX}_sorted.bam.bai"
fi
if [ -f "${PEAKS_FILE}.sorted" ] && [ "$PEAKS_FILE" != "$1" ]; then
    rm "${PEAKS_FILE}.sorted"
fi

echo "Analysis complete. Output files:"
echo "- ${OUTPUT_PREFIX}_peak_counts.txt"
echo "- ${OUTPUT_PREFIX}_viral_coverage.txt"
echo "- ${OUTPUT_PREFIX}_viral_depth.txt" 