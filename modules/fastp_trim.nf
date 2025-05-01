/*
 * define the TRIM process that trims the reads
 * given the read pairs
 */
process TRIM{
	fair true
	container 'staphb/fastp:0.24.0'
	debug true
	tag "fastp on $sample_id"
	publishDir	"${projectDir}/trimmed", mode: 'copy'

	maxForks 1

	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("*{1,2}.fastq.gz"), emit: trimmed_reads
 
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