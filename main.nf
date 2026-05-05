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
params.starindex = "$projectDir/star_index"
params.trimmeddir = "$projectDir/trimmed"
params.aligneddir = "$projectDir/aligned"
params.gtf = "$projectDir/gencode.v47.primary_assembly.basic.annotation.gtf"


/*
 * define the RNASEQ workflow that performs the RNAseq analysis
 * given the project directory
 */

workflow RNASEQ {
	refgenome = file("${projectDir}/*.genome.fa")

	if (params.single_end) {
		Channel
			.fromPath("${projectDir}/data/*.fastq.gz", checkIfExists: true)
			.map { file -> tuple(file.simpleName, [file]) }
			.set { read_pairs_ch }
	} else {
		Channel
			.fromFilePairs("${projectDir}/data/*{1,2}*.fastq.gz", checkIfExists: true)
			.set { read_pairs_ch }
	}
	read_pairs_ch.view()

	STAR_INDEX(refgenome)
	STAR_INDEX.out.view()


	// FASTQC(read_pairs_ch)
	// FASTQC.out.view()

	TRIM(read_pairs_ch)
	TRIM.out.view()

	ALIGN(STAR_INDEX.out, TRIM.out)
	ALIGN.out.view()

	FEATURECOUNT(params.gtf, ALIGN.out.bam.collect())
	FEATURECOUNT.out.view()

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
