/*
 * define the MULTIQC process that aggregates QC reports
 * Inputs: fastp JSON reports, STAR Log.final.out files, featureCounts summary
 */
process MULTIQC {
	container 'multiqc/multiqc:pdf-v1.34'
	publishDir "${params.outdir}/multiqc_report", mode: 'copy'

	input:
	path '*'

	output:
	path 'multiqc_report'

	script:
	"""
	multiqc . -o multiqc_report
	"""
}
