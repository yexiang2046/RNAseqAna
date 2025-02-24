process SPLICEQ {
	debug true
	container 'xiang2019/spliceq:v1.0.0'
	
	input:
	path bamfile
	path annotation

	output:
	path "*.tsv"

	script:
	"""
	bash bin/spliceq.sh $bamfile $annotation
	"""
	
}
