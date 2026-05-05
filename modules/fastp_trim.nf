/*
 * TRIM process - fastp quality control and trimming
 * Best practices:
 * - Automatic adapter detection and trimming
 * - Quality filtering (Q > 10 for alignment compatibility)
 * - PolyG tail trimming (for NovaSeq/NextSeq)
 * - Length filtering (min 36bp after trimming)
 * - 3' quality trimming with sliding window
 */
process TRIM {
	tag "fastp on $sample_id"
	publishDir "${params.outdir}/trimmed", mode: 'copy'

	container 'staphb/fastp:0.24.0'

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_*.fastq.gz"), emit: trimmed_reads
	path "${sample_id}_fastp.json", emit: json
	path "${sample_id}_fastp.html", emit: html

	script:
	def in2_arg   = params.single_end ? "" : "--in2 ${reads[1]}"
	def out2_arg  = params.single_end ? "" : "--out2 ${sample_id}_2.fastq.gz"
	def pe_args   = params.single_end ? "" : "--detect_adapter_for_pe --correction"
	"""
	fastp \
		--in1 ${reads[0]} \
		${in2_arg} \
		--out1 ${sample_id}_1.fastq.gz \
		${out2_arg} \
		--json ${sample_id}_fastp.json \
		--html ${sample_id}_fastp.html \
		--thread ${task.cpus} \
		${pe_args} \
		--qualified_quality_phred 10 \
		--unqualified_percent_limit 40 \
		--n_base_limit 5 \
		--length_required 36 \
		--cut_tail \
		--cut_tail_window_size 4 \
		--cut_tail_mean_quality 20 \
		--trim_poly_g \
		--trim_poly_x \
		--overrepresentation_analysis
	"""
}
