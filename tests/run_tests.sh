#!/bin/bash

# Exit on error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

echo "Starting pipeline tests..."

# Create test directories
mkdir -p test_results

# Run the test workflow
nextflow run tests/test_pipeline.nf \
    -profile test \
    --outdir test_results \
    -resume

# Check if test results exist
if [ ! -d "test_results" ]; then
    echo -e "${RED}Test failed: No results directory created${NC}"
    exit 1
fi

# Validate test outputs
echo "Validating test outputs..."

# Check for BAM preprocessing results
if [ ! -f "test_results/processed_bam/test.bam" ]; then
    echo -e "${RED}Test failed: Processed BAM file not found${NC}"
    exit 1
fi

# Check for peak calling results
if [ ! -f "test_results/peaks/test_peaks.bed" ]; then
    echo -e "${RED}Test failed: Peak file not found${NC}"
    exit 1
fi

# Check for feature annotation results
if [ ! -f "test_results/peak_annotations/test_feature_overlap_summary.csv" ]; then
    echo -e "${RED}Test failed: Feature annotation summary not found${NC}"
    exit 1
fi

# Check for repeat annotation results
if [ ! -f "test_results/repeat_annotations/test_rmsk_summary.csv" ]; then
    echo -e "${RED}Test failed: Repeat annotation summary not found${NC}"
    exit 1
fi

# Check for report generation
if [ ! -f "test_results/reports/summary_report.html" ]; then
    echo -e "${RED}Test failed: Summary report not found${NC}"
    exit 1
fi

echo -e "${GREEN}All tests passed successfully!${NC}" 