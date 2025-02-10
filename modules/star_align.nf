/*
 * define the STAR_INDEX process that creates a binary index
 * given the genome file
 */
process STAR_INDEX {
	debug true
	publishDir	"${params.output}/star_index", mode: 'copy'
	memory '32 GB'
	cpus 4
	containerOptions '-shm-size 32gb'

	input:
	path refgenome

	output:
	path "star_index"

	script:
	"""
	mkdir -p star_index
	STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${refgenome} --runThreadN ${params.cpus}
	"""
}


process ALIGN{
	fair true
	debug true
	tag "STAR on $sample_id"
	memory '32 GB'
	cpus 4
	containerOptions '-shm-size 32gb'

	publishDir "${params.output}/aligned", mode: 'copy'

	maxForks 1

	input:
	path star_index
	tuple   val(sample_id), path(reads)

	output:
	path "*Aligned.sortedByCoord.out.bam"

	script:
	"""
	STAR --genomeDir ${star_index} --readFilesIn ${reads[0]} ${reads[1]} \
    	--readFilesCommand zcat --runThreadN ${params.cpus} --genomeLoad NoSharedMemory      \
    	--outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    	--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
    	--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    	--outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    	--outSAMtype BAM SortedByCoordinate --sjdbScore 1     \
	--limitBAMsortRAM 1000409716 --outFileNamePrefix ${sample_id}
	"""
}
