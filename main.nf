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


params.metadata = "$projectDir/metadata.txt"


/*
 * define the RNASEQ workflow that performs the RNAseq analysis
 * given the project directory
 */



workflow RNASEQ {
	refgenome = file("${projectDir}/*.genome.fa")	

	
	Channel
	   	.fromFilePairs("${projectDir}/data/*{1,2}*.fastq.gz", checkIfExists: true)
	   	.set { read_pairs_ch }
	read_pairs_ch.view()

	STAR_INDEX(refgenome)
	STAR_INDEX.out.view()


	// FASTQC(read_pairs_ch)
	// FASTQC.out.view()
	
	TRIM(read_pairs_ch)
	TRIM.out.view()

	ALIGN(STAR_INDEX.out, TRIM.out)
	ALIGN.out.view()

	// Run Piranha peak calling on aligned BAM files
	PIRANHA_PEAK_CALLING(ALIGN.out)
	PIRANHA_PEAK_CALLING.out.peaks.view()

	FEATURECOUNT(params.gtf, ALIGN.out.collect())

	FEATURECOUNT.out.view()

	// Define input channels
    Channel.fromPath('counts.txt').set { counts_ch }
    Channel.fromPath('metadata.txt').set { metadata_ch }

    // Run DE analysis
    DE_ANALYSIS(FEATURECOUNT.out, params.metadata)
    // DE_ANALYSIS(counts_ch, metadata_ch)

    // Run Functional Analysis
    // FUNCTIONAL_ANALYSIS(de_results_ch)
    // FUNCTIONAL_ANALYSIS(DE_ANALYSIS.out)

	// emit: FASTQC.out | concat(TRIM.out) | concat(ALIGN.out) | collect	
	// emit: FASTQC.out | concat(TRIM.out) | concat(ALIGN.out) | collect	
}

workflow  {
	RNASEQ()
}

// Process to run Piranha peak calling
process PIRANHA_PEAK_CALLING {
    publishDir "${params.outdir}/piranha_output/\${bam.simpleName}", mode: 'copy'
    
    input:
    path bam
    
    output:
    tuple val(bam.simpleName), path("*_peaks.bed"), emit: peaks
    path "*_converted.bed", emit: converted_bed
    
    script:
    def piranha_opts = params.piranha_params ? params.piranha_params : ''
    """
    # Convert BAM to BED
    bedtools bamtobed -i $bam > ${bam.simpleName}_converted.bed
    
    # Sort the BED file
    sort -k1,1 -k2,2n ${bam.simpleName}_converted.bed > ${bam.simpleName}_sorted.bed
    
    # Run Piranha on the sorted BED file
    piranha $piranha_opts ${bam.simpleName}_sorted.bed > ${bam.simpleName}_peaks.bed
    
    # Clean up intermediate files
    rm ${bam.simpleName}_sorted.bed
    """
}
