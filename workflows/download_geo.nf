#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DOWNLOAD_GEO } from '../modules/download_geo'

workflow {
    // Input parameters
    params.geo_accession = null
    params.outdir = 'results'
    
    // Validate required parameters
    if (!params.geo_accession) {
        error "Please provide a GEO accession ID using --geo_accession"
    }
    
    // Create output directory
    mkdir(params.outdir)
    
    // Run the download process
    DOWNLOAD_GEO(
        params.geo_accession,
        params.outdir
    )
    
    // Print summary
    DOWNLOAD_GEO.out.fastq_files.collect().subscribe { fastq_files ->
        log.info """
        ==============================================
        Download completed successfully!
        
        Output files:
        - Metadata: ${params.outdir}/${params.geo_accession}_metadata.txt
        - Sample info: ${params.outdir}/sample_info.txt
        - edgeR metadata: ${params.outdir}/metadata.csv
        - Contrasts: ${params.outdir}/contrasts.csv
        - Raw reads: ${params.outdir}/raw_reads/
        
        Number of FASTQ files: ${fastq_files.size()}
        ==============================================
        """
    }
} 