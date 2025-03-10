#!/usr/bin/env nextflow

/*
 * Salmon quantification workflow
 */

// Parameters
params.reads = "$projectDir/data/*_{1,2}.fastq.gz"
params.transcriptome = "$projectDir/reference/transcriptome.fa"
params.genome = "$projectDir/reference/genome.fa"
params.outdir = "results"
params.metadata = "$projectDir/metadata.csv"
params.contrasts = "$projectDir/contrasts.csv"

// Log info
log.info """\
    SALMON QUANTIFICATION WORKFLOW
    =============================
    reads        : ${params.reads}
    transcriptome: ${params.transcriptome}
    genome      : ${params.genome}
    outdir      : ${params.outdir}
    metadata    : ${params.metadata}
    contrasts   : ${params.contrasts}
    """
    .stripIndent()

// Import modules
include { SALMON_INDEX } from '../modules/salmon_index'
include { SALMON_QUANT } from '../modules/salmon_quant'
include { FASTQC } from '../modules/fastqc'
include { MULTIQC } from '../modules/multiqc'
include { DE_ANALYSIS } from '../modules/de_analysis'

workflow {
    // Create channels
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    
    transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
    
    // Quality control
    FASTQC(read_pairs_ch)
    
    // Create Salmon index
    SALMON_INDEX(transcriptome_ch, genome_ch)
    
    // Run Salmon quantification
    SALMON_QUANT(read_pairs_ch, SALMON_INDEX.out.index)
    
    // MultiQC report
    MULTIQC(
        FASTQC.out.zip_files.collect()
    )
    
    // Differential expression analysis
    DE_ANALYSIS(
        SALMON_QUANT.out.quant_results.collect(),
        file(params.metadata),
        file(params.contrasts)
    )
} 