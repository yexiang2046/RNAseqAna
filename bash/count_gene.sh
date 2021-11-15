#!/bin/bash

featureCounts -p -a wgEncodeGencodeBasicV34.gtf -o unique_counts_gene.txt -g gene_id -B -P -C -T 28 aligned/731-JK-1_S19Aligned.sortedByCoord.out.bam aligned/731-JK-2_S20Aligned.sortedByCoord.out.bam aligned/731-JK-3_S21Aligned.sortedByCoord.out.bam aligned/731-JK-4_S22Aligned.sortedByCoord.out.bam aligned/731-JK-5_S23Aligned.sortedByCoord.out.bam aligned/731-JK-6_S24Aligned.sortedByCoord.out.bam aligned/731-JK-7_S25Aligned.sortedByCoord.out.bam aligned/731-JK-8_S26Aligned.sortedByCoord.out.bam
