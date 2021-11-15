#!/bin/bash





export PATH=$PATH:~/src/STAR-2.7.3a/bin/Linux_x86_64_static/
GRCh38=$1
idx=$2
ncpus=$3

STAR --runMode genomeGenerate --genomeDir ${idx} --genomeFastaFiles ${GRCh38} --runThreadN ${ncpus}
