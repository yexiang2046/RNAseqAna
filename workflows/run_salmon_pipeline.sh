#!/bin/bash

# Default parameters
THREADS=12
LIBTYPE="A"  # Automatic detection of library type

# Help message
usage() {
    echo "Usage: $0 [-h] -r READ_DIR -t TRANSCRIPTOME -g GENOME -o OUTPUT_DIR -m METADATA [-c CONTRASTS] [-p THREADS]"
    echo "  -r READ_DIR      Directory containing paired-end reads (*_{1,2}.fastq.gz)"
    echo "  -t TRANSCRIPTOME Path to transcriptome FASTA file"
    echo "  -g GENOME        Path to genome FASTA file (for decoy-aware index)"
    echo "  -o OUTPUT_DIR    Output directory"
    echo "  -m METADATA      Metadata CSV file"
    echo "  -c CONTRASTS     Contrasts CSV file (optional)"
    echo "  -p THREADS       Number of threads (default: 12)"
    echo "  -h              Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hr:t:g:o:m:c:p:" opt; do
    case $opt in
        h) usage ;;
        r) READ_DIR="$OPTARG" ;;
        t) TRANSCRIPTOME="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) METADATA="$OPTARG" ;;
        c) CONTRASTS="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$READ_DIR" ] || [ -z "$TRANSCRIPTOME" ] || [ -z "$GENOME" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$METADATA" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"/{fastqc,salmon_index,salmon_quant,multiqc,de_analysis}

# Run FastQC
echo "Running FastQC..."
fastqc -t "$THREADS" "$READ_DIR"/*_{1,2}.fastq.gz -o "$OUTPUT_DIR/fastqc"

# Create Salmon index
echo "Creating Salmon index..."
# Generate decoys file
grep "^>" "$GENOME" | cut -d " " -f 1 | sed 's/>//g' > "$OUTPUT_DIR/salmon_index/decoys.txt"
# Concatenate transcriptome and genome
cat "$TRANSCRIPTOME" "$GENOME" > "$OUTPUT_DIR/salmon_index/gentrome.fa"

# Build index
salmon index \
    -t "$OUTPUT_DIR/salmon_index/gentrome.fa" \
    -i "$OUTPUT_DIR/salmon_index/salmon_index" \
    -d "$OUTPUT_DIR/salmon_index/decoys.txt" \
    -p "$THREADS"

# Run Salmon quantification
echo "Running Salmon quantification..."
for r1 in "$READ_DIR"/*_1.fastq.gz; do
    r2="${r1/_1./_2.}"
    sample_name=$(basename "$r1" _1.fastq.gz)
    
    echo "Processing sample: $sample_name"
    salmon quant \
        -i "$OUTPUT_DIR/salmon_index/salmon_index" \
        -l "$LIBTYPE" \
        -1 "$r1" \
        -2 "$r2" \
        --validateMappings \
        --gcBias \
        --seqBias \
        -o "$OUTPUT_DIR/salmon_quant/$sample_name" \
        -p "$THREADS"
done

# Run MultiQC
echo "Running MultiQC..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR/multiqc"

# Run differential expression analysis
echo "Running differential expression analysis..."
if [ -n "$CONTRASTS" ]; then
    Rscript bin/edger.r "$OUTPUT_DIR/salmon_quant" "$METADATA" "$CONTRASTS"
else
    Rscript bin/edger.r "$OUTPUT_DIR/salmon_quant" "$METADATA"
fi

echo "Pipeline completed successfully!"
echo "Results can be found in: $OUTPUT_DIR" 