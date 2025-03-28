#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DOWNLOAD_GEO } from '../modules/download_geo'

workflow {
    // Input parameters
    params.id_list = null
    params.outdir = 'results'
    
    // Validate required parameters
    if (!params.id_list) {
        error "Please provide an input file containing SRA IDs using --id_list"
    }
    
    // Create output directory
    def outdir = file(params.outdir)
    outdir.mkdirs()
    
    // Run the download process
    DOWNLOAD_GEO(
        file(params.id_list),
        outdir
    )
    
    // Print summary
    DOWNLOAD_GEO.out.fastq_files.collect().subscribe { fastq_files ->
        log.info """
        ==============================================
        Download completed successfully!
        
        Output files:
        - Raw reads: ${outdir}/raw_reads/
        
        Number of FASTQ files: ${fastq_files.size()}
        ==============================================
        """
    }
} 