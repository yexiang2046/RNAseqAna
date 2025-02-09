process MULTIQC {
	debug true
	publishDir "${params.output}/multiqc_out", mode:'copy'

	input:
	path '*'

	output:
	path 'multiqc_report'

	script:
	"""
	multiqc .
	"""
}