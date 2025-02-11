/*
 * define the FEATURECOUNT process that performs feature counting
 * given the GTF file and the BAM file
 */
process FEATURECOUNT {
	debug true
	container 'biocontainers/subread:v1.6.3dfsg-1-deb_cv1'
	publishDir "${projectDir}", mode:'copy'

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