
process FEATURECOUNT {
	debug true
	publishDir "${params.output}", mode:'copy'
	containerOptions '-shm-size 10gb'

	input:
	val    gtf
    path    bamfile

    output:
    path    "*.txt"

	script:
	"""
	featureCounts -T 14 -p -t exon -g gene_id -F GTF -a ${projectDir}/${gtf} -o counts.txt ${bamfile}
	"""
}
