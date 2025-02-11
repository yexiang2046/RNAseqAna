/*
 * define the TRIM process that trims the reads
 * given the read pairs
 */
process TRIM{
	fair true
	container 'biocontainers/fastp:v0.20.1_cv1'
	debug true
	tag "fastp on $sample_id"
	publishDir	"${params.trimmeddir}", mode: 'copy'

	maxForks 3

	input:
	tuple	val(sample_id), path(reads)

	output:
	tuple val("$sample_id"),path("*{1,2}.fastp.fastq.gz")
 
	script:
	"""
	fastp -w 16 -l 20 -i ${reads[0]} -I ${reads[1]} -o ${sample_id}1.fastp.fastq.gz -O ${sample_id}2.fastp.fastq.gz
	"""
}