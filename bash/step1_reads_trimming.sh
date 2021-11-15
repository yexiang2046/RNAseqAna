#!/bin/bash


threads=$1
align_idx=$2
out_dir=$3
trimmer_dir=$4
alignment_dir=$5
pairend=$6
trimmed_dir=$7
fastq_file1=$8
fastq_file2=$9

export PATH=$PATH:$trimmer_dir

if [[ $pairend == 'Yes' ]]; then
  echo "Trim paired reads: $pairend"
  trim_galore --illumina --gzip --paired -q 20 --fastqc -j $threads $fastq_file1 $fastq_file2 -o $trimmed_dir
else
  echo "Trim single reads: No"
  trim_galore --illumina --gzip -q 20 --fastqc -j $threads $fastq_file1 -o $trimmed_dir
fi
