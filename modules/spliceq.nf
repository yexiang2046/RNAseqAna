process SPLICEQ {
	debug true
	container 'xiang2019/spliceq:v1.0.0'
	publishDir "${params.projectDir/spliceq_out}", mode:'copy'
	
	input:
	path(bamfile)
	path annotation
	val(sample_id)

	output:
	path "*.tsv"

	script:
	"""
	bash spliceq.sh $bamfile $annotation $sample_id
	"""
	
}
