#!/bin/bash


export PATH=$PATH:/home/xiang/src/FastQC

FASTQC_OUT=/home/xiang/DNA_seq/BJAB_LANA_RNAseq/fastqc_out
THREADS=8
FASTQ_DIR=/home/xiang/DNA_seq/BJAB_LANA_RNAseq/OG

mkdir $FASTQC_OUT
fastqc -o $FASTQC_OUT -t $THREADS -f fastq $FASTQ_DIR/*
