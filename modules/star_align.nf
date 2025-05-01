/*
 * define the STAR_INDEX process that creates a binary index
 * given the genome file
 */
process STAR_INDEX {
	container 'alexdobin/star:2.6.1d'
	publishDir "${projectDir}", mode: 'copy'

	input:
	path refgenome

	output:
	path "star_index", emit: star_index

	script:
	"""
	mkdir -p star_index
	STAR --runMode genomeGenerate \
		--genomeDir star_index \
		--genomeFastaFiles ${refgenome} \
		--runThreadN ${params.cpus}
	"""
}

/*
 * define the ALIGN process that aligns the reads
 * given the trimmed read pairs and the STAR index
 */
process ALIGN {
	container 'alexdobin/star:2.6.1d'
	tag "STAR on $sample_id"
	publishDir "${projectDir}/aligned", mode: 'copy'
	maxForks 1
	input:
	path star_index
	tuple   val(sample_id), path(reads)

	output:
	path("*Aligned.sortedByCoord.out.bam"), emit: bam

	script:
	"""
	STAR --genomeDir ${star_index} \
		--readFilesIn ${reads[0]} ${reads[1]} \
		--readFilesCommand zcat \
		--runThreadN ${params.cpus} \
		--genomeLoad NoSharedMemory \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--outFilterType BySJout \
		--outSAMattributes NH HI AS NM MD \
		--outSAMtype BAM SortedByCoordinate \
		--sjdbScore 1 \
		--limitBAMsortRAM ${params.ram} \
		--outFileNamePrefix ${sample_id}
	"""
}