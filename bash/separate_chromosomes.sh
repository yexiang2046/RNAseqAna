#!/bin/bash

# Usage: bash separate_chromosomes.sh your_fasta_file.fasta

fasta_file=$1

# Check if the input file is provided
if [ -z "$fasta_file" ]; then
    echo "Usage: bash separate_chromosomes.sh your_fasta_file.fasta"
    exit 1
fi

# Use awk to split the file into separate chromosome files
awk '/^>/ {OUT=substr($0,2) ".fasta"} OUT {print > OUT}' $fasta_file

echo "Chromosomes have been separated into individual FASTA files."
