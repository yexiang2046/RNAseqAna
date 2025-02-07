process MULTIQC {
	debug true
	publishDir "${params.projectDir}", mode:'copy'

	input:
	path '*'

	output:
	path 'multiqc_report'

	script:
	"""
	multiqc .
	"""
}