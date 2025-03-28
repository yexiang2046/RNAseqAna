process DOWNLOAD_GEO {
    container "ncbi/sra-tools:3.1.0"
    
    input:
        path geo_accession
        path output_dir
    
    output:
        path "${output_dir}/${geo_accession}_metadata.txt", emit: metadata
        path "${output_dir}/sample_info.txt", emit: sample_info
        path "${output_dir}/metadata.csv", emit: edger_metadata
        path "${output_dir}/contrasts.csv", emit: contrasts
        path "${output_dir}/raw_reads/*.fastq", emit: fastq_files
    
    script:
    """
    #!/usr/bin/env bash

    
    # Run the script
    bin/download_geo_data.sh -g ${geo_accession} -o ${output_dir}
    """
} 