#!/usr/bin/env nextflow

params.cpus = 12 
params.ram = 60000000000 /* ~60GB */

params.starindex = "$projectDir/star_index"
params.trimmeddir = "$projectDir/trimmed"

params.aligneddir = "$projectDir/aligned"
/* params.refgenome = "$projectDir/GRCh38.primary_assembly.genome.fa" */
params.gtf = "$projectDir/gencode.v47.primary_assembly.annotation.gtf"

params.refgenomelink = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz"
params.refgtflink = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz"


process FASTQC {
	debug true
	tag "FASTQC on $sample_id"

	container 'nfcore/rnaseq:1.4.2'

	input:
	tuple val(sample_id), path(reads)

	output:
	path "fastqc_${sample_id}_logs"

	script:
	"""
	mkdir "fastqc_${sample_id}_logs"
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""
}
	
/*
 * define the STAR_INDEX process that creates a binary index
 * given the genome file
 */
process STAR_INDEX {
	debug true
	publishDir	"${params.starindex}", mode: 'copy'

	input:
	path refgenome
	path starindex

	output:
	path $starindex

	script:
	"""
	mkdir ${params.starindex}
	STAR --runMode genomeGenerate --genomeDir ${starindex} --genomeFastaFiles ${refgenome} --runThreadN ${params.cpus}
	"""
}

process TRIM{
	debug true
	tag "fastp on $sample_id"
	publishDir	"${params.trimmeddir}", mode: 'copy'

	input:
	tuple	val(sample_id), path(reads)

	output:
	tuple	val(sample_id), path("${sample_id}1.fastp.fastq.gz"), path("${sample_id}*2.fastp.fastq.gz")

	script:
	"""
	fastp -w 16 -l 20 -i ${reads[0]} -I ${reads[1]} -o ${sample_id}1.fastp.fastq.gz -O ${sample_id}2.fastp.fastq.gz
	"""
}



process ALIGN{
	debug true
	tag "STAR on $sample_id"

	publishDir "${params.aligneddir}", mode: 'copy'

	input:
	path star_index
	tuple	val(sample_id), path(read1), path(read2) 

	output:
	path	"${sample_id}.bam"

	script:
	"""
	STAR --genomeDir ${star_index} --readFilesIn ${read1} ${read2} \
    		--readFilesCommand zcat --runThreadN ${params.cpus} --genomeLoad NoSharedMemory      \
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

process DESEQ2_QC {
	debug true
	publishDir "${projectDir}/deseq2_qc_out"
	input:
	path val(count_file) /* output from commandline featureCounts */

	output:
	path "deseq2_qc_out"

	script:
	"""
		deseq2_qc.r -i ${count_file} -f 7 -r "Aligned.sortedByCoord.out.bam" -o "deseq2_qc_out" -c ${params.cpus}
	"""
}

workflow {
	Channel
		.fromPath("${projectDir}/*.genome.fa")
		.set{refgenome_ch}

	
	Channel
	   		.fromFilePairs("${projectDir}/data/*{1,2}_001.fastq.gz", checkIfExists: true)
	   		.set { read_pairs_ch }
	read_pairs_ch.view()

	index_ch = STAR_INDEX( refgenome_ch, params.starindex )
	index_ch.view()


	fastqc_ch = FASTQC(read_pairs_ch)
	fastqc_ch.view()
	
	reads_ch = TRIM(read_pairs_ch)
	reads_ch.view()

	align_ch = ALIGN(index_ch, reads_ch)
	align_ch.view()

    MULTIQC(align_ch.mix(fastqc_ch).collect())

	

	bamfile_ch = align_ch.filter("*.bam")
	bamfile_ch.view()

	FEATURECOUNT(bamfile_ch)

	
}
