#!/bin/bash


# Configuration variables
STAR_BIN_DIR=~/src/STAR-2.7.3a/bin/Linux_x86_64_static
THREADS=8
LIMIT_RAM=50000000000

# Set Bash options for error handling and command execution
set -Eeuo pipefail


# Function to display help message
display_help() {
    echo "Usage: $0 FILE_PATH READ1_PAT READ2_PAT STAR_IDX OUT_PATH"
    echo "Arguments:"
    echo "  FILE_PATH    Path to the directory containing fastq.gz files"
    echo "  READ1_PAT    Pattern for Read 1 fastq.gz files"
    echo "  READ2_PAT    Pattern for Read 2 fastq.gz files"
    echo "  STAR_IDX     Path to the STAR index directory"
    echo "  OUT_PATH     Output directory path"
    echo "  GENOME	 GENOME(fasta) to align to"
}

# Function to log errors
log_error() {
    echo "Error: $1" >&2
}

# Trap errors and log them
trap 'log_error "Script error in line $LINENO"; exit 1' ERR


# Check for the help option
if [[ "$1" == "--help" ]]; then
    display_help
    exit 0
fi

export PATH=$PATH:${STAR_BIN_DIR}

# assign input variables
FILE_PATH="$1"
READ1_PAT="$2"
READ2_PAT="$3"
STAR_IDX="$4"
OUT_PATH="$5"
GENOME="$6"

# Check if the input file paths exist
if [[ ! -d "$FILE_PATH" ]]; then
    log_error "Input directory $FILE_PATH does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUT_PATH"

# generate index
mkdir -p ${STAR_IDX}
STAR --runMode genomeGenerate --genomeDir ${STAR_IDX} --genomeFastaFiles ${GENOME} --runThreadN 12

if [[ ! -f "$STAR_IDX/Genome" ]]; then
    log_error "STAR index directory $STAR_IDX is invalid or incomplete."
    exit 1
fi



samples_id=$(find "$FILE_PATH" -maxdepth 1 -type f -name '*.fastq.gz' -printf '%P\n' | awk -F'.' '{print $1}' | sort | uniq)
# samples_id=$(ls ${FILE_PATH} | cut -d '.' -f 1 | sort | uniq)
printf '%s\n' "$samples_id" | while IFS= read -r f; do
    echo "Read1: $FILE_PATH/$f${READ1_PAT}"
    echo "Read2: $FILE_PATH/$f${READ2_PAT}"
    STAR --genomeDir "$STAR_IDX" --readFilesIn "$FILE_PATH/$f${READ1_PAT}" "$FILE_PATH/$f${READ2_PAT}" \
    --readFilesCommand zcat --runThreadN ${THREADS} --genomeLoad NoSharedMemory      \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate --sjdbScore 1     \
    --limitBAMsortRAM ${LIMIT_RAM} --outFileNamePrefix "$OUT_PATH/$f"
done
chmod -R a+wr ${OUT_PATH}
