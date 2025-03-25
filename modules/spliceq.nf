process SPLICEQ {
	debug true
	publishDir "${projectDir}/spliceq_out", mode:'copy'
	
	input:
	tuple val(sample_id), path(bamfile)
	path annotation
	

	output:
	path "*.tsv"

	script:
	"""
	spliceq.sh $bamfile $annotation $sample_id
	"""
	
}
