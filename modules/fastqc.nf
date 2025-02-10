
process FASTQC {
	debug true
	tag "FASTQC on $sample_id"
	containerOptions '-shm-size 60gb'
	memory '64 GB'
	publishDir "${params.output}/fastqc_out", mode:'copy'
	maxForks 4

	input:
	tuple val(sample_id), path(reads)

	output:
	path "fastqc_*_logs"

	script:
	"""
	mkdir "fastqc_${sample_id}_logs"
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""
}
	
