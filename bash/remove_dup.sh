#!/bin/bash

export PATH=$PATH:~/src/STAR-2.7.3a/bin/Linux_x86_64_static
for bam in $(ls aligned | grep out.bam)
do
	STAR --inputBAMfile ${bam} \
	--bamRemoveDuplicatesType UniqueIdentical \
	--runMode inputAlignmentsFromBAM \
	--bamRemoveDuplicatesMate2basesN 15 \
	--outFileNamePrefix dupremoved.${bam} \
	--limitBAMsortRAM 30000000000
done 
