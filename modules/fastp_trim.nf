process TRIM{
	fair true
	debug true
	tag "fastp on $sample_id"
	memory '16 GB'  
    	cpus 4 
	containerOptions '-shm-size 10gb'
	publishDir	"${params.output}/trimmed", mode: 'copy'

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
