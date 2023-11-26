#!/bin/bash

set -Eeuo pipefail


# folder for fastq files
FILE_PATH=$1
ALIGNMENT_PATH=$2
# can be mouse, human or custom chromsomes provided
CHROMOSOMES=$3


# mouse chromosomes
MOUSE_CHROMOSOM="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"
HUMAN_CHROMOSOM="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# assign CHROMSOMES
if [ ${CHROMOSOMES} == "mouse" ]; then
    CHROMOSOMES_KEPT=${MOUSE_CHROMOSOM}
elif [ ${CHROMOSOMES} == "human" ]; then
    CHROMOSOMES_KEPT=${HUMAN_CHROMOSOM}
else 
    CHROMOSOMES_KEPT=${CHROMOSOMES}
fi



samples_id=$(find "$FILE_PATH" -maxdepth 1 -type f -name '*.fastq.gz' -printf '%P\n' | awk -F'_' '{print $1}' | sort | uniq)

for f in ${samples_id}
do
    docker run -v $(pwd):$(pwd) -v ${FILE_PATH}:${FILE_PATH} -v ${ALIGNMENT_PATH}:${ALIGNMENT_PATH} -w $(pwd) xiang2019/samtools:v1.9 samtools view -@ 8 -q 255 -b -h -o ${ALIGNMENT_PATH}/${f}.unique.bam ${ALIGNMENT_PATH}/${f}.bamAligned.sortedByCoord.out.bam ${CHROMOSOMES_KEPT}
done


