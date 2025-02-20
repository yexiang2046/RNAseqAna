process FASTQC {
	debug true
	tag "FASTQC on $sample_id"

	container 'xiang2019/rnaseq_cmd:v1.0.0'
	publishDir "${projectDir}/fastqc_logs", mode:'copy'
	maxForks 1

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