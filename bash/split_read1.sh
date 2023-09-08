#!/bin/bash

set -Eeuo pipefail


# folder for fastq files
FILE_PATH=$1
ALIGNMENT_PATH=$2




samples_id=$(ls ${FILE_PATH} | awk -v FS='\t' -v ORS='\n' 'match($0, /^.*_001.fastq.gz$/) {print substr($0, RSTART, RLENGTH-24)}' | awk '!seen[$0]++')


for f in ${samples_id}
do
    docker run -v $(pwd):$(pwd) -v ${FILE_PATH}:${FILE_PATH} -v ${ALIGNMENT_PATH}:${ALIGNMENT_PATH} -w $(pwd) xiang2019/samtools:v1.9 samtools view -@ 8 -f 64 -b -h ${ALIGNMENT_PATH}/${f}.unique.bam > ${ALIGNMENT_PATH}/${f}.unique.read1.bam
done