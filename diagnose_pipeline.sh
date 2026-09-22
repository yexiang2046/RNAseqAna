#!/bin/bash

echo "=== RNAseq Pipeline Diagnostics ==="
echo ""

# Check for Nextflow log
echo "1. Checking for Nextflow logs..."
if [ -f .nextflow.log ]; then
    echo "✓ Found .nextflow.log"
    echo ""
    echo "Last 50 lines of log:"
    tail -50 .nextflow.log
else
    echo "✗ No .nextflow.log found - are you in the directory where you ran Nextflow?"
fi

echo ""
echo "2. Checking reference genome files..."
GENOME_FILES=$(ls *.genome.fa 2>/dev/null)
if [ -n "$GENOME_FILES" ]; then
    echo "✓ Found reference genome: $GENOME_FILES"
    ls -lh *.genome.fa
else
    echo "✗ No *.genome.fa file found in current directory"
fi

echo ""
echo "3. Checking input FASTQ files..."
if [ -d "data" ]; then
    echo "✓ data/ directory exists"
    FASTQ_COUNT=$(ls data/*.fastq.gz 2>/dev/null | wc -l)
    echo "   Found $FASTQ_COUNT FASTQ files:"
    ls -lh data/*.fastq.gz 2>/dev/null || echo "   No .fastq.gz files in data/"
else
    echo "✗ data/ directory not found"
fi

echo ""
echo "4. Checking work directory..."
if [ -d "work" ]; then
    echo "✓ work/ directory exists"
    # Find STAR_INDEX work directory
    STAR_INDEX_WORK=$(find work -type f -name ".command.sh" -exec grep -l "STAR.*genomeGenerate" {} \; | head -1 | xargs dirname 2>/dev/null)
    if [ -n "$STAR_INDEX_WORK" ]; then
        echo "   STAR_INDEX work directory: $STAR_INDEX_WORK"
        echo ""
        echo "   Checking STAR_INDEX outputs:"
        ls -la "$STAR_INDEX_WORK/"
        echo ""
        echo "   STAR_INDEX stdout (.command.out):"
        cat "$STAR_INDEX_WORK/.command.out" 2>/dev/null | tail -20
        echo ""
        echo "   STAR_INDEX stderr (.command.err):"
        cat "$STAR_INDEX_WORK/.command.err" 2>/dev/null | tail -20
        echo ""
        echo "   STAR_INDEX exit status:"
        cat "$STAR_INDEX_WORK/.exitcode" 2>/dev/null
    else
        echo "   Could not find STAR_INDEX work directory"
    fi
else
    echo "✗ work/ directory not found"
fi

echo ""
echo "5. Checking results/outputs..."
if [ -d "results" ]; then
    echo "✓ results/ directory exists"
    echo "   Contents:"
    find results -type f 2>/dev/null | head -20
else
    echo "✗ No results/ directory found"
fi

echo ""
echo "6. Memory and resource check..."
echo "   Available memory:"
free -h
echo ""
echo "   CPU info:"
nproc
echo ""
echo "   Docker status:"
docker ps -a --filter "name=nxf-" --format "table {{.ID}}\t{{.Image}}\t{{.Status}}" 2>/dev/null || echo "   Docker not accessible or no Nextflow containers"

echo ""
echo "=== End Diagnostics ==="
