#!/bin/bash

# add fastqc to path
export PATH=$PATH:/home/xiang/src/FastQC

# help info
if [ $# -lt 3 ]; then
    echo "Usage: $0 <fastqc_out_dir> <threads> <fastq_file>"
    echo "Example: $0 fastqc_out 4 sample.fastq"
    echo "Note: fastqc_out_dir will be created if not exists"
    exit 1
fi

# parameters
FASTQC_OUT=$1
THREADS=$2
FASTQ_FILE=$3


mkdir ${FASTQC_OUT}
fastqc -o ${FASTQC_OUT} -t ${THREADS} -f fastq ${FASTQ_FILE}
