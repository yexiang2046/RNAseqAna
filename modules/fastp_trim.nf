/*
 * define the TRIM process that trims the reads
 * given the read pairs
 */
process TRIM {
	tag "fastp on $sample_id"
	publishDir "${projectDir}/trimmed", mode: 'copy'
	
	container 'staphb/fastp:0.24.0'
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("*_{1,2}.fastq.gz"), emit: trimmed_reads
 
	script:
	"""
	fastp \
		--in1 ${reads[0]} \
		--in2 ${reads[1]} \
		--out1 ${sample_id}_1.fastq.gz \
		--out2 ${sample_id}_2.fastq.gz \
		--thread ${task.cpus} \
		--detect_adapter_for_pe \
		--cut_mean_quality 20 \
		--cut_window_size 4
	"""
}