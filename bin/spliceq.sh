#!/bin/bash

# get the input arguments
bamfile=$1
annotation=$2
sample_id=$3
# run the spliceq analysis
# output default tsv file
SPLICE-q.py -b $bamfile -g $annotation -o ${sample_id}.csv

