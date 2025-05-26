/*
 * define the FEATURECOUNT process that performs feature counting
 * given the GTF file and the BAM file
 */
process FEATURECOUNT {
	debug true
	container 'xiang2019/rnaseq_cmd:v1.0.0'
	publishDir "${projectDir}/feature_counts", mode:'copy'

	input:
	path    gtf
    path    bamfile

    output:
    path    "*.txt"

	script:
	"""
	featureCounts -T 14 -p -t exon -g gene_id -F GTF -a ${gtf} -o counts.txt ${bamfile}
	"""
}