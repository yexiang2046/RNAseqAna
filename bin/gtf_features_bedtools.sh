#!/bin/bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i INPUT_GTF -o OUTPUT_DIR"
    echo "Extract genomic features from GTF file using bedtools"
    echo ""
    echo "Arguments:"
    echo "  -i INPUT_GTF     Input GTF file"
    echo "  -o OUTPUT_DIR    Output directory for BED files"
    echo "  -h              Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:o:" opt; do
    case $opt in
        h) usage ;;
        i) INPUT="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$INPUT" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input file exists
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file does not exist: $INPUT"
    exit 1
fi

# Check if bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Processing features using bedtools..."

# First, extract protein-coding transcripts with transcript IDs
echo "Extracting protein-coding transcripts..."
awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ && $3=="transcript" {
    split($9, attrs, "; ")
    is_coding = 0
    transcript_id = ""
    gene_name = ""
    
    for (i in attrs) {
        if (attrs[i] ~ /gene_type "protein_coding"/ || attrs[i] ~ /transcript_type "protein_coding"/) {
            is_coding = 1
        }
        if (attrs[i] ~ /transcript_id/) {
            split(attrs[i], tid, " ")
            transcript_id = tid[2]
        }
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gname, " ")
            gene_name = gname[2]
        }
    }
    
    if (is_coding && transcript_id != "") {
        gsub(/"/, "", transcript_id)
        gsub(/"/, "", gene_name)
        print $1"\t"($4-1)"\t"$5"\t"transcript_id"\t"$6"\t"$7"\t"gene_name
    }
}' "$INPUT" > "$OUTPUT_DIR/temp_protein_coding_transcripts.bed"

# Extract CDS regions with transcript IDs and gene names (with duplicate removal)
echo "Extracting CDS regions..."
awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ && $3=="CDS" {
    split($9, attrs, "; ")
    transcript_id = ""
    gene_name = ""
    gene_id = ""
    
    for (i in attrs) {
        if (attrs[i] ~ /transcript_id/) {
            split(attrs[i], tid, " ")
            transcript_id = tid[2]
        }
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gname, " ")
            gene_name = gname[2]
        }
        if (attrs[i] ~ /gene_id/) {
            split(attrs[i], gid, " ")
            gene_id = gid[2]
        }
    }
    
    if (transcript_id != "") {
        gsub(/"/, "", transcript_id)
        gsub(/"/, "", gene_name)
        gsub(/"/, "", gene_id)
        print $1"\t"($4-1)"\t"$5"\t"transcript_id"\t"$6"\t"$7"\t"gene_name"\t"gene_id
    }
}' "$INPUT" | sort -k1,1 -k2,2n -k3,3n -k4,4 -k7,7 | uniq > "$OUTPUT_DIR/cds.bed"

# Extract exons with gene names (with duplicate removal)
echo "Extracting exons..."
awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ && $3=="exon" {
    split($9, attrs, "; ")
    transcript_id = ""
    gene_name = ""
    gene_id = ""
    
    for (i in attrs) {
        if (attrs[i] ~ /transcript_id/) {
            split(attrs[i], tid, " ")
            transcript_id = tid[2]
        }
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gname, " ")
            gene_name = gname[2]
        }
        if (attrs[i] ~ /gene_id/) {
            split(attrs[i], gid, " ")
            gene_id = gid[2]
        }
    }
    
    if (transcript_id != "") {
        gsub(/"/, "", transcript_id)
        gsub(/"/, "", gene_name)
        gsub(/"/, "", gene_id)
        print $1"\t"($4-1)"\t"$5"\t"transcript_id"\t"$6"\t"$7"\t"gene_name"\t"gene_id
    }
}' "$INPUT" | sort -k1,1 -k2,2n -k3,3n -k4,4 -k7,7 | uniq > "$OUTPUT_DIR/exons.bed"

# Extract genes with gene names (with duplicate removal)
echo "Extracting genes..."
awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ && $3=="gene" {
    split($9, attrs, "; ")
    gene_id = ""
    gene_name = ""
    gene_type = ""
    
    for (i in attrs) {
        if (attrs[i] ~ /gene_id/) {
            split(attrs[i], gid, " ")
            gene_id = gid[2]
        }
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gname, " ")
            gene_name = gname[2]
        }
        if (attrs[i] ~ /gene_type/) {
            split(attrs[i], gtype, " ")
            gene_type = gtype[2]
        }
    }
    
    if (gene_id != "") {
        gsub(/"/, "", gene_id)
        gsub(/"/, "", gene_name)
        gsub(/"/, "", gene_type)
        print $1"\t"($4-1)"\t"$5"\t"gene_id"\t"$6"\t"$7"\t"gene_name"\t"gene_type
    }
}' "$INPUT" | sort -k1,1 -k2,2n -k3,3n -k4,4 -k7,7 | uniq > "$OUTPUT_DIR/genes.bed"

# Process UTRs with gene names
echo "Processing UTRs..."
awk -F'\t' 'BEGIN {OFS="\t"}
    # First pass: store CDS boundaries for each transcript
    NR==FNR {
        transcript_id = $4
        gene_name = $7
        gene_id = $8
        if (!(transcript_id in cds_start) || $2 < cds_start[transcript_id]) {
            cds_start[transcript_id] = $2
            cds_chr[transcript_id] = $1
            cds_strand[transcript_id] = $6
            cds_gene_name[transcript_id] = gene_name
            cds_gene_id[transcript_id] = gene_id
        }
        if (!(transcript_id in cds_end) || $3 > cds_end[transcript_id]) {
            cds_end[transcript_id] = $3
        }
        next
    }
    # Second pass: process transcripts and output UTRs
    {
        transcript_id = $4
        if (transcript_id in cds_start) {
            if ($6 == "+") {
                # 5UTR
                if ($2 < cds_start[transcript_id]) {
                    print $1, $2, cds_start[transcript_id], \
                          transcript_id"_5UTR", ".", $6, \
                          cds_gene_name[transcript_id], \
                          cds_gene_id[transcript_id] > "'"$OUTPUT_DIR"'/five_prime_utr.bed"
                }
                # 3UTR
                if ($3 > cds_end[transcript_id]) {
                    print $1, cds_end[transcript_id], $3, \
                          transcript_id"_3UTR", ".", $6, \
                          cds_gene_name[transcript_id], \
                          cds_gene_id[transcript_id] > "'"$OUTPUT_DIR"'/three_prime_utr.bed"
                }
            } else if ($6 == "-") {
                # 5UTR
                if ($3 > cds_end[transcript_id]) {
                    print $1, cds_end[transcript_id], $3, \
                          transcript_id"_5UTR", ".", $6, \
                          cds_gene_name[transcript_id], \
                          cds_gene_id[transcript_id] > "'"$OUTPUT_DIR"'/five_prime_utr.bed"
                }
                # 3UTR
                if ($2 < cds_start[transcript_id]) {
                    print $1, $2, cds_start[transcript_id], \
                          transcript_id"_3UTR", ".", $6, \
                          cds_gene_name[transcript_id], \
                          cds_gene_id[transcript_id] > "'"$OUTPUT_DIR"'/three_prime_utr.bed"
                }
            }
        }
    }' "$OUTPUT_DIR/cds.bed" "$OUTPUT_DIR/temp_protein_coding_transcripts.bed"

# Sort and remove duplicates from UTR files
sort -u -k1,1 -k2,2n -k3,3n -k4,4 -k7,7 -o "$OUTPUT_DIR/five_prime_utr.bed" "$OUTPUT_DIR/five_prime_utr.bed" 
sort -u -k1,1 -k2,2n -k3,3n -k4,4 -k7,7 -o "$OUTPUT_DIR/three_prime_utr.bed" "$OUTPUT_DIR/three_prime_utr.bed"

# Clean up temporary files
rm "$OUTPUT_DIR/temp_"*.bed

# Print statistics
echo -e "\nFeature counts:"
echo "Protein-coding transcripts: $(wc -l < "$OUTPUT_DIR/temp_protein_coding_transcripts.bed") transcripts"
echo "CDS regions: $(wc -l < "$OUTPUT_DIR/cds.bed") regions"
echo "5' UTRs: $(wc -l < "$OUTPUT_DIR/five_prime_utr.bed") regions"
echo "3' UTRs: $(wc -l < "$OUTPUT_DIR/three_prime_utr.bed") regions"

echo -e "\nFiles created in $OUTPUT_DIR:"
ls -l "$OUTPUT_DIR" 