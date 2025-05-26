#!/bin/bash

# Check if folder path is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <folder_path>"
    exit 1
fi

FOLDER_PATH="$1"

# Check if folder exists
if [ ! -d "$FOLDER_PATH" ]; then
    echo "Error: $FOLDER_PATH is not a valid directory"
    exit 1
fi

# Process each pair of fastq files
for f1 in "$FOLDER_PATH"/*_R1*.fastq.gz; do
    # Get corresponding R2 file
    f2=${f1/_R1/_R2}
    
    # Check if both files exist
    if [ ! -f "$f2" ]; then
        echo "Warning: Missing pair for $f1"
        continue
    fi

    # Get base name for output files
    base=$(basename "$f1" | sed 's/_R1.*/.random1M/')
    out1="${base}_R1.fastq.gz"
    out2="${base}_R2.fastq.gz"

    echo "Processing pair: $f1 and $f2"
    
    # Process the pair
    paste <(zcat "$f1") <(zcat "$f2") | \
        awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | \
        shuf | \
        head -n 1000000 | \
        sed 's/\t\t/\n/g' | \
        awk -v out1=>(gzip > "$out1") -v out2=>(gzip > "$out2") \
            '{print $1 > out1; print $2 > out2}'

    echo "Created: $out1 and $out2"
done

echo "Done! All pairs processed"

