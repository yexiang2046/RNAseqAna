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
			.fromPath("${params.data_dir}/*.fastq.gz", checkIfExists: true)
			.map { file -> tuple(file.simpleName, [file]) }
			.set { read_pairs_ch }
	} else {
		Channel
			.fromFilePairs("${params.data_dir}/*{1,2}*.fastq.gz", checkIfExists: true)
			.set { read_pairs_ch }
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
