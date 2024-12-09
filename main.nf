#!/usr/bin/env nextflow

params.cpus = 8
params.ram = 60000000000 /* ~60GB */



params.trimeddir = "$projectDir/trimmed"
params.aligneddir = "$projectDir/aligned"
params.refgenome = "$projectDir/GRCh38.primary_assembly.genome.fa"
params.gtf = "$projectDir/gencode.v47.primary_assembly.annotation.gtf"

params.refgenomelink = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz"
params.refgtflink = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz"


process FASTQC {
	debug true
	tag "FASTQC on $sample_id"

	input:
	tuple val(sample_id), path(reads)

	output:
	path "fastqc_${sample_id}_logs"

	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""
}
	
/*
 * define the STAR_INDEX process that creates a binary index
 * given the genome file
 */
process STAR_INDEX {
	debug true

	input:
	path refgenome

	output:
	path 'star_index'

	script:
	"""
	STAR --runMode genomeGenerate --genomeDir ${$params.starindex} --genomeFastaFiles ${refgenome} --runThreadN ${params.cpus}
	"""
}

process TRIM{
	debug true
	tag "fastp on $sample_id"

	input:
	tuple	val(sample_id), path(reads)

	output:
	path	"*.fastp.fastq.gz"
	script:
	"""
	fastp -w 16 -l 20 -i ${reads[0]} -I ${reads[1]} -o ${params.trimmeddir}/${sample_id} -O ${params.trimmeddir}/${sample_id}
	"""
}



process ALIGN{
	debug true
	tag "STAR on $sample_id"

	publishDir params.aligneddir, mode: 'copy'

	input:
	path star_index
	tuple	val(sample_id), path(reads) 

	output:
	path	"${sample_id}.bam"

	script:
	"""
	STAR --genomeDir ${star_index} --readFilesIn ${reads[0]} ${reads[1]} \
    		--readFilesCommand zcat --runThreadN ${cpus} --genomeLoad NoSharedMemory      \
    		--outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    		--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
    		--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    		--outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    		--outSAMtype BAM SortedByCoordinate --sjdbScore 1     \
    		--limitBAMsortRAM ${params.ram} --outFileNamePrefix ${params.aligneddir}/${sample_id}
	"""
}

process MULTIQC {
	debug true
	publishDir "${params.projectDir}", mode:'copy'

	input:
	path '*'

	output:
	path 'multiqc_report.html'

	script:
	"""
	multiqc .
	"""
}

process FEATURECOUNT {
	debug true
	publishDir "$params.projectDir}/featureCounts", mode: 'copy'

	input:
	tuple	val(sample_id), path(bamfile)


	script:
	"""
	featureCounts -T 14 -p -t exon -g gene_id -F GTF -a ${params.gtf} -o ${sample_id}.txt ${bamfile}
	"""
}

workflow {
	index_ch = STAR_INDEX(params.refgenome)
	index_ch.view()

	Channel
	   		.fromFilePairs("${projectDir}/data/*{1,2}_001.fastq.gz", checkIfExists: true)
	   		.set { read_pairs_ch }
	read_pairs_ch.view()

	fastqc_ch = FASTQC(read_pairs_ch)
	fastqc_ch.view()
	
	trim_read_pairs_ch = TRIM(read_pairs_ch)
	trim_read_pairs_ch.view()

	align_ch = ALIGN(index_ch, trim_read_pairs_ch)
	align_ch.view()

    MULTIQC(align_ch.mix(fastqc_ch).collect())

	

	bamfile_ch = align_ch.filter("*.bam")
	bamfile_ch.view()

	FEATURECOUNT(bamfile_ch)

	
}
