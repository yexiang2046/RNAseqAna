/*
 * define the TRIM process that trims the reads
 * given the read pairs
 */
process TRIM{
	fair true
	container 'xiang2019/rnaseq_cmd:v1.0.0'
	debug true
	tag "fastp on $sample_id"
	publishDir	"${projectDir}/trimmed", mode: 'copy'

	maxForks 1

	input:
	tuple	val(sample_id), path(reads)

	output:
	tuple val("$sample_id"),path("*{1,2}.fastp.fastq.gz")
 
	script:
	"""
	fastp -w 16 -l 20 -i ${reads[0]} -I ${reads[1]} -o ${sample_id}1.fastp.fastq.gz -O ${sample_id}2.fastp.fastq.gz
	"""
}