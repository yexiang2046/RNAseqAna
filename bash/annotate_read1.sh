#!/bin/bash

set -Eeuo pipefail


# folder for fastq files
FILE_PATH=$1
ALIGNMENT_PATH=$2




samples_id=$(ls ${FILE_PATH} | awk -v FS='\t' -v ORS='\n' 'match($0, /^.*_001.fastq.gz$/) {print substr($0, RSTART, RLENGTH-24)}' | awk '!seen[$0]++')


for f in ${samples_id}
do
    docker run -v $(pwd):$(pwd) -v ${FILE_PATH}:${FILE_PATH} -v ${ALIGNMENT_PATH}:${ALIGNMENT_PATH} -w $(pwd) quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0 bedtools intersect -wao -split -bed -a ${ALIGNMENT_PATH}/${f}.unique.read1.bam -b hg38_rmsk_gencodev43.gtf | awk -e '$13 ~ /chr/ {print $0}' - > ${ALIGNMENT_PATH}/${f}.anno.bed
    sort -k 4 ${ALIGNMENT_PATH}/${f}.anno.bed > ${ALIGNMENT_PATH}/${f}.anno.sorted.bed
    awk 'BEGIN { FS = "\t"; OFS = " " };{gsub(/gene_id|transcript_type|transcript_id/, "")}; /hg38_rmsk/ {split($21, TR, ";"); READ1[$4] = READ1[$4] TR[1]}; /rtracklayer/ {split($21, TR, ";"); READ1[$4] = READ1[$4] $15 TR[3]};  END {for (NAME in READ1) print NAME, READ1[NAME]}' ${ALIGNMENT_PATH}/${f}.anno.sorted.bed > ${ALIGNMENT_PATH}/${f}.txt
    awk '{ while(++i<=NF) printf (!a[$i]++) ? $i FS : ""; i=split("",a); print "" }' ${ALIGNMENT_PATH}/${f}.txt > ${ALIGNMENT_PATH}/${f}.dedup.txt
    awk -f bash/rank_ann.awk ${ALIGNMENT_PATH}/${f}.dedup.txt| cut -d ' ' -f 3 | sort | uniq -c > ${ALIGNMENT_PATH}/${f}.tx.type.txt
done