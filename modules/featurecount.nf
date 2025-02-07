
process FEATURECOUNT {
	debug true
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