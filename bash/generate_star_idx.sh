#!/bin/bash


## genome=resources/GRCh38.primary_GQ994935.fa


export PATH=$PATH:~/src/STAR-2.7.3a/bin/Linux_x86_64_static/
genome=$1
idx=$2
ncpus=$3

STAR --runMode genomeGenerate --genomeDir ${idx} --genomeFastaFiles ${genome} --runThreadN ${ncpus}
