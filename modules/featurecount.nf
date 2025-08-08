/*
 * define the FEATURECOUNT process that performs feature counting
 * given the GTF file and the BAM file
 */
process FEATURECOUNT {
	debug true
	cpus = 8
	container 'xiang2019/rnaseq_cmd:v1.0.0'
	publishDir "${params.outdir}/feature_counts", mode:'copy'

	input:
	path gtf
	path bamfile

	output:
	path "counts.txt", emit: counts
	path "counts.txt.summary", emit: summary

	script:
	"""
	featureCounts -T ${task.cpus} -p -t exon -g gene_id -F GTF -a ${gtf} -o counts.txt ${bamfile}
	"""
}
