#!/bin/bash

export PATH=$PATH:~/src/STAR-2.7.3a/bin/Linux_x86_64_static/
ram_GB=50
echo "using $1 threads:"
ncpus=$1
echo "The star index folder is $2:"
star_index_folder=$2
echo "The current directory is: $3"
current_dir=$3
echo "The fastq files used are: $4 and $5 as pairs"
FASTQ1=$4
FASTQ2=$5

sample=$(basename $FASTQ1)

# outputfiles
ALIGNMENT_DIR=$6
echo "the alignment files are in $ALIGNMENT_DIR"

#mkdir ${star_index_folder}
#STAR --runMode genomeGenerate --genomeDir ${star_index_folder} --genomeFastaFiles ./GRCh38.primary_GQ994935.fa --runThreadN ${ncpus}
## align with ENCODE STANDARD option
if [[ "$FASTQ2" == "No" ]]
then
  STAR --genomeDir $star_index_folder --readFilesIn ${FASTQ1} \
    --readFilesCommand zcat --runThreadN $ncpus --genomeLoad NoSharedMemory      \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate --sjdbScore 1     \
    --limitBAMsortRAM ${ram_GB}000000000 --outFileNamePrefix ${ALIGNMENT_DIR}/${sample}

else

  STAR --genomeDir $star_index_folder --readFilesIn ${FASTQ1} ${FASTQ2} \
    --readFilesCommand zcat --runThreadN $ncpus --genomeLoad NoSharedMemory      \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate --sjdbScore 1     \
    --limitBAMsortRAM ${ram_GB}000000000 --outFileNamePrefix ${ALIGNMENT_DIR}/${sample}

fi
