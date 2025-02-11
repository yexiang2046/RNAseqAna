#!/usr/bin/env nextflow

/*
 * import the modules
 */
include {
	modules.fastqc
	modules.fastp_trim
	modules.star_align
	modules.featurecount
	modules.de_analysis
	modules.multiqc
}

/*
 * define the parameters
 */
params.cpus = 12 
params.ram = 60000000000 /* ~60GB */


params.projectDir = "/home/xiang/Projects/RNAseqAna"
params.starindex = "$projectDir/star_index"
params.trimmeddir = "$projectDir/trimmed"
params.aligneddir = "$projectDir/aligned"



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

	ALIGN(STAR_INDEX.out.collect(), TRIM.out)
	ALIGN.out.view()

	FEATURECOUNT(params.gtf, ALIGN.out.collect())

	FEATURECOUNT.out.view()

	// Define input channels
    Channel.fromPath('counts.txt').set { counts_ch }
    Channel.fromPath('metadata.txt').set { metadata_ch }
    Channel.fromPath('de_results').set { output_ch }

    // Run DE analysis
    DE_ANALYSIS(counts_ch, metadata_ch, output_ch)

	emit: FASTQC.out | concat(TRIM.out) | concat(ALIGN.out) | collect	
}

workflow  {
	RNASEQ()
}
