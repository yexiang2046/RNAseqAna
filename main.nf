#!/usr/bin/env nextflow

/*
 * import the modules
 */
include { FASTQC } from './modules/fastqc.nf'
include { STAR_INDEX; ALIGN } from './modules/star_align.nf'
include { TRIM } from './modules/fastp_trim.nf'
include { FEATURECOUNT } from './modules/featurecount.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { DE_ANALYSIS } from './modules/de_analysis.nf'
include { FUNCTIONAL_ANALYSIS } from './modules/functional_analysis.nf'

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


params.metadata = "$projectDir/metadata.txt"


/*
 * define the RNASEQ workflow that performs the RNAseq analysis
 * given the project directory
 */



workflow RNASEQ {
	refgenome = file("${projectDir}/*.genome.fa")	

	
	Channel
	   	.fromFilePairs("${projectDir}/data/*{1,2}_001.fastq.gz", checkIfExists: true)
	   	.set { read_pairs_ch }
	read_pairs_ch.view()

	STAR_INDEX(refgenome)
	STAR_INDEX.out.view()


	FASTQC(read_pairs_ch)
	FASTQC.out.view()
	
	TRIM(read_pairs_ch)
	TRIM.out.view()

	ALIGN(STAR_INDEX.out, TRIM.out)
	ALIGN.out.view()

	FEATURECOUNT(params.gtf, ALIGN.out.collect())

	FEATURECOUNT.out.view()


    // Run DE analysis
    DE_ANALYSIS(FEATURECOUNT.out, params.metadata)

    // Run Functional Analysis
    // FUNCTIONAL_ANALYSIS(de_results_ch)

	// emit: FASTQC.out | concat(TRIM.out) | concat(ALIGN.out) | collect	
}

workflow  {
	RNASEQ()
}
