#!/bin/bash


export PATH=$PATH:/home/xiang/src/FastQC

FASTQC_OUT=$1
THREADS=$2
FASTQ_FILE=$3


mkdir ${FASTQC_OUT}
fastqc -o ${FASTQC_OUT} -t ${THREADS} -f fastq ${FASTQ_FILE}
