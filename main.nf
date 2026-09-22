#!/usr/bin/env nextflow

/*
 * import the modules
 */
include { FASTQC } from './modules/fastqc.nf'
include { STAR_INDEX; ALIGN } from './modules/star_align.nf'
include { TRIM } from './modules/fastp_trim.nf'
include { FEATURECOUNT } from './modules/featurecount.nf'
include { MULTIQC } from './modules/multiqc.nf'


/*
 * define the parameters
 */
params.cpus = 12
params.ram = 60000000000 /* ~60GB */


params.projectDir = "/home/xiang/Projects/RNAseqAna"
params.data_dir   = "${params.projectDir}/data"
params.gtf        = "${params.projectDir}/gencode.v47.primary_assembly.basic.annotation.gtf"


/*
 * define the RNASEQ workflow that performs the RNAseq analysis
 * given the project directory
 */

workflow RNASEQ {
	refgenome = file("${projectDir}/*.genome.fa")

	if (params.single_end) {
		Channel
			.fromPath("${params.data_dir}/*.{fastq,fq}.gz", checkIfExists: true)
			.map { file -> tuple(file.simpleName, [file]) }
			.set { read_pairs_ch }
	} else {
		// If user provides a custom pattern, use it
		if (params.read_pattern) {
			Channel
				.fromFilePairs("${params.data_dir}/${params.read_pattern}", checkIfExists: true, size: 2)
				.set { read_pairs_ch }
		} else {
			// Auto-detect common paired-end naming patterns
			// Try patterns in order of specificity
			def patterns = [
				['*_R{1,2}_*.{fastq,fq}.gz', 'Illumina style with lane: sample_R1_001.fastq.gz'],
				['*_R{1,2}.{fastq,fq}.gz', 'Underscore R1/R2: sample_R1.fastq.gz'],
				['*_{1,2}.{fastq,fq}.gz', 'Underscore 1/2: sample_1.fastq.gz'],
				['*.R{1,2}.{fastq,fq}.gz', 'Dot R1/R2: sample.R1.fastq.gz'],
				['*R{1,2}.{fastq,fq}.gz', 'No separator R1/R2: sampleR1.fastq.gz'],
				['*{1,2}.{fastq,fq}.gz', 'Just numbers: sample1.fastq.gz'],
				['*_{1,2}_*.{fastq,fq}.gz', 'Underscore with lane: sample_1_001.fastq.gz']
			]
			
			// Test which pattern finds files
			def matched_pattern = null
			for (pattern_info in patterns) {
				def pattern = pattern_info[0]
				def test_glob = "${params.data_dir}/${pattern}"
				def test_files = file(test_glob)
				
				if (test_files && (test_files instanceof List ? test_files.size() > 0 : true)) {
					matched_pattern = pattern
					log.info "✓ Detected paired-end files using pattern: ${pattern_info[1]}"
					log.info "  Pattern: ${pattern}"
					break
				}
			}
			
			if (matched_pattern) {
				Channel
					.fromFilePairs("${params.data_dir}/${matched_pattern}", checkIfExists: true, size: 2)
					.set { read_pairs_ch }
			} else {
				error """
No paired-end FASTQ files found in ${params.data_dir}/
Tried common patterns:
  - sample_R1.fastq.gz / sample_R2.fastq.gz
  - sample_1.fastq.gz / sample_2.fastq.gz  
  - sample.R1.fastq.gz / sample.R2.fastq.gz
  
Please ensure your files follow one of these naming conventions,
or specify a custom pattern with --read_pattern

Examples:
  --read_pattern '*_R{1,2}.fastq.gz'
  --read_pattern '*_{1,2}.fq.gz'
"""
			}
		}
	}

	if (params.star_index) {
		star_index_ch = Channel.value(file(params.star_index, checkIfExists: true))
	} else {
		STAR_INDEX(refgenome)
		star_index_ch = STAR_INDEX.out.star_index
	}

	// FASTQC(read_pairs_ch)

	TRIM(read_pairs_ch)

	ALIGN(star_index_ch, TRIM.out.trimmed_reads)

	FEATURECOUNT(file(params.gtf), ALIGN.out.bam.collect())

	MULTIQC(
		TRIM.out.json
			.mix(ALIGN.out.log)
			.mix(FEATURECOUNT.out.summary)
			.collect()
	)
}

workflow  {
	RNASEQ()
}
