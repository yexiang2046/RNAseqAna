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
include { PIRANHA_PEAK_CALLING } from './modules/piranha_peak_calling.nf'


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
params.genome = "$projectDir/GRCh38.primary_assembly_KSHV.genome.fa"
params.metadata = "$projectDir/metadata.txt"


/*
 * define the RNASEQ workflow that performs the RNAseq analysis
 * given the project directory
 */



workflow RNASEQ {
	refgenome = file(params.genome)	

	// Create channel for read pairs
	Channel
	   	.fromFilePairs("${projectDir}/data/*{1,2}*.fastq.gz", checkIfExists: true)
	   	.set { read_pairs_ch }

	// Run STAR indexing
	STAR_INDEX(refgenome)
	
	// Run trimming
	TRIM(read_pairs_ch)
	
	// Run alignment
	ALIGN(STAR_INDEX.out, TRIM.out.trimmed_reads)
	
	// Run feature counting
	FEATURECOUNT(params.gtf, ALIGN.out.bam.collect())

	// Run DE analysis
	DE_ANALYSIS(FEATURECOUNT.out, params.metadata)
}

workflow {
	RNASEQ()
}


