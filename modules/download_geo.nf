process DOWNLOAD_GEO {
    container "xiang2019/sratools:v3.2.1"
    
    input:
        path id_list
        path output_dir
    
    output:
        path "${output_dir}/raw_reads/*.fastq", emit: fastq_files
    
    script:
    """
    # Run the script
    ${baseDir}/../bin/download_SRA_data.sh -i ${id_list} -o ${output_dir}
    """
} 