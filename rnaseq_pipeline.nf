#!/usr/bin/env nextflow

/*
 * import the modules
 */
include { FASTQC } from './modules/fastqc.nf'
include { STAR_INDEX; ALIGN } from './modules/star_align.nf'
include { TRIM } from './modules/fastp_trim.nf'
include { FEATURECOUNT } from './modules/featurecount.nf'
include { DE_ANALYSIS } from './modules/de_analysis.nf'
include { FUNCTIONAL_ANALYSIS } from './modules/functional_analysis.nf'


/*
 * define the parameters
 */
params.cpus = 12 
params.ram = 60000000000 /* ~60GB */

params.readspath = "data"
params.outdir = "results"
params.starindex = "$projectDir/star_index"
params.trimmeddir = "$projectDir/trimmed"
params.aligneddir = "$projectDir/aligned"
params.gtf = "$projectDir/gencode.v47.primary_assembly.basic.annotation.gtf"
params.genome = "$projectDir/GRCh38.primary_assembly_KSHV.genome.fa"
params.metadata = "$projectDir/metadata.txt"
params.species = "human"



/*
 * define the RNASEQ workflow that performs the RNAseq analysis
 * given the project directory
 */



workflow RNASEQ {
	refgenome = file(params.genome)	

	// Validate required files exist
	if (!file(params.genome).exists()) {
		error "Genome file not found: ${params.genome}"
	}
	if (!file(params.gtf).exists()) {
		error "GTF file not found: ${params.gtf}"
	}
	if (!file(params.metadata).exists()) {
		error "Metadata file not found: ${params.metadata}"
	}

	// Create channel for read pairs - more flexible pattern
        Channel
                .fromFilePairs("${params.readspath}/*{_R1,_R2}_*.fastq.gz", checkIfExists: true)
                .ifEmpty { error "No read pairs found in ${params.readspath} matching pattern *{_R1,_R2}_*.fastq.gz" }
                .set { read_pairs_ch }

        // Debug output for read pairs
        println "DEBUG: Read pairs channel contents:"
        read_pairs_ch.view()

	
	// Debug output for read pairs
	println "DEBUG: Read pairs channel contents:"
	read_pairs_ch.view()

	// Run STAR indexing
	STAR_INDEX(refgenome)
	
	// Debug output for STAR index
	println "DEBUG: STAR index output:"
	STAR_INDEX.out.star_index.view()
	
	// Run trimming
	TRIM(read_pairs_ch)
	
	// Debug output for trimmed reads
	println "DEBUG: Trimmed reads output:"
	TRIM.out.trimmed_reads.view()
	
	// Run alignment
	ALIGN(STAR_INDEX.out.star_index, TRIM.out.trimmed_reads)
	
	// Debug output for aligned BAMs
	println "DEBUG: Aligned BAMs output:"
	ALIGN.out.bam.view()
	
	// Run feature counting
	FEATURECOUNT(file(params.gtf), ALIGN.out.bam.collect())
	
	// Run deg analysis
	DE_ANALYSIS(FEATURECOUNT.out.counts, file(params.metadata), file(params.gtf), params.species)

	// Run functional analysis on DEG results (if enabled)
	if (params.run_functional_analysis) {
		// Transform DEG files to tuple with comparison name
		DE_ANALYSIS.out.deg_files
			.flatten()
			.map { file ->
				def comparison_name = file.name.replaceAll(/^DEG_(.+)\.csv$/, '$1')
				tuple(comparison_name, file)
			}
			.set { deg_files_ch }

		FUNCTIONAL_ANALYSIS(deg_files_ch)
	}

}

workflow {
	RNASEQ()
}


