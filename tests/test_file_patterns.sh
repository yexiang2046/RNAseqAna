#!/bin/bash
# Test script to verify FASTQ file pattern detection

set -e

echo "=== Testing FASTQ File Pattern Detection ==="
echo ""

# Create a temporary test directory
TEST_DIR=$(mktemp -d)
echo "Test directory: $TEST_DIR"
echo ""

# Function to test a specific pattern
test_pattern() {
    local pattern_name=$1
    local file1=$2
    local file2=$3
    
    echo "Testing pattern: $pattern_name"
    echo "  Files: $file1, $file2"
    
    # Create test data directory
    mkdir -p "$TEST_DIR/data"
    
    # Create dummy FASTQ files
    echo "@read1" | gzip > "$TEST_DIR/data/$file1"
    echo "@read2" | gzip > "$TEST_DIR/data/$file2"
    
    # Create a minimal test workflow
    cat > "$TEST_DIR/test.nf" <<'EOF'
#!/usr/bin/env nextflow

params.data_dir = "data"
params.read_pattern = null

workflow {
    if (params.read_pattern) {
        Channel
            .fromFilePairs("${params.data_dir}/${params.read_pattern}", checkIfExists: true, size: 2)
            .view { sample_id, files -> "Detected: ${sample_id} with files: ${files}" }
    } else {
        def patterns = [
            ['*_R{1,2}_*.{fastq,fq}.gz', 'Illumina with lane'],
            ['*_R{1,2}.{fastq,fq}.gz', 'Underscore R1/R2'],
            ['*_{1,2}.{fastq,fq}.gz', 'Underscore 1/2'],
            ['*.R{1,2}.{fastq,fq}.gz', 'Dot R1/R2'],
            ['*R{1,2}.{fastq,fq}.gz', 'No separator R1/R2'],
            ['*{1,2}.{fastq,fq}.gz', 'Just numbers'],
            ['*_{1,2}_*.{fastq,fq}.gz', 'Underscore with lane']
        ]
        
        def matched_pattern = null
        for (pattern_info in patterns) {
            def pattern = pattern_info[0]
            def test_glob = "${params.data_dir}/${pattern}"
            def test_files = file(test_glob)
            
            if (test_files && (test_files instanceof List ? test_files.size() > 0 : true)) {
                matched_pattern = pattern
                log.info "✓ Detected using: ${pattern_info[1]}"
                break
            }
        }
        
        if (matched_pattern) {
            Channel
                .fromFilePairs("${params.data_dir}/${matched_pattern}", checkIfExists: true, size: 2)
                .view { sample_id, files -> "Detected: ${sample_id} with files: ${files}" }
        } else {
            error "No files found"
        }
    }
}
EOF
    
    # Run the test
    cd "$TEST_DIR"
    if nextflow run test.nf 2>&1 | grep -q "Detected:"; then
        echo "  ✓ Pattern detected successfully"
    else
        echo "  ✗ Pattern detection failed"
    fi
    
    # Cleanup for next test
    rm -rf "$TEST_DIR/data" "$TEST_DIR/test.nf" "$TEST_DIR/.nextflow"* "$TEST_DIR/work"
    echo ""
}

# Test common patterns
test_pattern "Illumina R1/R2 with lane" "sample_R1_001.fastq.gz" "sample_R2_001.fastq.gz"
test_pattern "Underscore R1/R2" "sample_R1.fastq.gz" "sample_R2.fastq.gz"
test_pattern "Underscore 1/2" "sample_1.fastq.gz" "sample_2.fastq.gz"
test_pattern "Dot R1/R2" "sample.R1.fastq.gz" "sample.R2.fastq.gz"
test_pattern "No separator R1/R2" "sampleR1.fastq.gz" "sampleR2.fastq.gz"
test_pattern "Just numbers" "sample1.fastq.gz" "sample2.fastq.gz"
test_pattern "fq.gz extension" "sample_R1.fq.gz" "sample_R2.fq.gz"

# Cleanup
rm -rf "$TEST_DIR"

echo ""
echo "=== All Pattern Tests Complete ==="
