#!/bin/bash


thread_num=$1
input_bam=$2

output_sense=$3
output_antisense=$4

bamCoverage -p $1 -b $2 -o $3 --binSize 10 --filterRNAstrand reverse --normalizeUsing CPM

bamCoverage -p $1 -b $2 -o $4 --binSize 10 --filterRNAstrand forward --normalizeUsing CPM
