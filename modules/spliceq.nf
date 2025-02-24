process SPLICEQ {
	debug true
	container 'xiang2019/spliceq:v1.0.0'
	publishDir "${projectDir}/spliceq_out", mode:'copy'
	
	input:
	tuple val(sample_id), path(bamfile)
	path annotation
	

	output:
	path "*.tsv"

	script:
	"""
	bash spliceq.sh $bamfile $annotation $sample_id
	"""
	
}
