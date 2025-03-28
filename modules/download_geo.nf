process DOWNLOAD_GEO {
    container "ncbi/sra-tools:3.1.0"
    
    input:
        path id_list
        path output_dir
    
    output:
        path "${output_dir}/raw_reads/*.fastq", emit: fastq_files
    
    script:
    """
    # Run the script
    bin/download_SRA_data.sh -i ${id_list} -o ${output_dir}
    """
} 