process MULTIQC {
	debug true
	container 'nfcore/rnaseq:1.4.2'
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