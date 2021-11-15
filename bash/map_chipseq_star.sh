#!/bin/bash


### This script mapped fastq to pre_index genome with STAR, tailed for ChIP-seq data (do not allow splicing of reads)
### Also generate bigwig file with MAPQ > 20, using bin window of 10bp

ncpus=30

idx=star_idx

# path for STAR
export PATH=$PATH:~/src/STAR-2.7.3a/bin/Linux_x86_64_static/



mkdir aligned
for sample in $(ls OG | grep fastq)
do
STAR --genomeDir ./${idx} --readFilesIn ./OG/${sample} --readFilesCommand zcat --runThreadN 30 --genomeLoad NoSharedMemory \
	--outFilterMultimapNmax 500 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
	--alignIntronMax 1 --alignMatesGapMax 1000000 --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterScoreMinOverLread 0.85 \
	--seedSearchStartLmax 30 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000 --outFileNamePrefix ./aligned/${sample}
done 



##### index bam files
for file in $(ls ./aligned | grep bam$)
do
	samtools index -@ 20 ./aligned/${file}
done

for file in $(ls ./aligned | grep out.bam)
do
	bamCoverage -b ./aligned/${file} -o ./aligned/${file}.bw -bs 10 --normalizeUsing CPM -p 30 --minMappingQuality 20 -of bigwig
done

