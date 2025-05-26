/*
 * define the MULTIQC process that performs multiqc analysis
 * given the project directory
 */
process MULTIQC {
	debug true
	container 'xiang2019/rnaseq_cmd:v1.0.0'
	publishDir "${params.projectDir}/multiqc_report", mode:'copy'

	input:
	path '*'

	output:
	path 'multiqc_report'

	script:
	"""
	multiqc .
	"""
}