#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.id_list = 'SRR_ID.txt'
params.output_dir = 'results'

process DOWNLOAD_GEO {
    container 'xiang2019/sratools:v3.2.1'
    
    publishDir 'results', mode: 'copy'

    input:
    path srr_file    // Input file containing SRR IDs

    output:
    path "results/raw_reads/*.fastq", emit: fastq_files

    script:
    """
    mkdir -p results/raw_reads
    
    
    # Run the script
    ${baseDir}/../bin/download_SRA_data.sh -i ${srr_file} -o ${params.output_dir}
    """
}

workflow {
    // Define the input channel
    srr_ids = Channel.fromPath(params.id_list)
    
    // Run the download process
    DOWNLOAD_GEO(srr_ids)
} 