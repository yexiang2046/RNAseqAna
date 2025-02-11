process FASTQC {
	debug true
	tag "FASTQC on $sample_id"

	container 'nfcore/rnaseq:1.4.2'

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