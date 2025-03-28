process DOWNLOAD_GEO {
    container "ncbi/sra-tools:3.1.0"
    
    input:
        val geo_accession
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
    
    # Download SRA metadata from Run Selector
    SRA_URL="https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${geo_accession}&output=file"
    curl -s "\$SRA_URL" > "${output_dir}/${geo_accession}_metadata.txt"
    
    # Extract SRA IDs and sample information
    awk -F',' '
        NR==1 {
            for(i=1; i<=NF; i++) {
                if(\$i=="Run") run_col=i
                if(\$i=="Sample Name") sample_col=i
                if(\$i=="condition") condition_col=i
                if(\$i=="cell_type") cell_type_col=i
                if(\$i=="tissue") tissue_col=i
                if(\$i=="time_point") time_point_col=i
            }
        }
        NR>1 {
            # Clean up the fields (remove quotes and escape characters)
            gsub(/"/, "", \$run_col)
            gsub(/"/, "", \$sample_col)
            gsub(/"/, "", \$condition_col)
            gsub(/"/, "", \$cell_type_col)
            gsub(/"/, "", \$tissue_col)
            gsub(/"/, "", \$time_point_col)
            gsub(/\\,/, "_", \$sample_col)
            gsub(/\\,/, "_", \$condition_col)
            gsub(/\\,/, "_", \$cell_type_col)
            gsub(/\\,/, "_", \$tissue_col)
            gsub(/\\,/, "_", \$time_point_col)
            
            # Combine condition, cell_type, tissue, and time_point for a more descriptive condition
            condition = \$condition_col
            if (\$cell_type_col != "") condition = condition "_" \$cell_type_col
            if (\$tissue_col != "") condition = condition "_" \$tissue_col
            if (\$time_point_col != "") condition = condition "_" \$time_point_col
            
            print \$run_col "\\t" \$sample_col "\\t" condition
        }
    ' "${output_dir}/${geo_accession}_metadata.txt" > "${output_dir}/sample_info.txt"
    
    # Create metadata.csv for edgeR
    echo "sample,condition" > "${output_dir}/metadata.csv"
    awk -F'\\t' '
        {
            print \$1 "," \$3
        }
    ' "${output_dir}/sample_info.txt" >> "${output_dir}/metadata.csv"
    
    # Create contrasts.csv for edgeR
    echo "name,treatment,control" > "${output_dir}/contrasts.csv"
    # Get unique conditions
    conditions=\$(awk -F',' 'NR>1 {print \$2}' "${output_dir}/metadata.csv" | sort -u)
    # Create contrasts between all pairs of conditions
    first_condition=\$(echo "\$conditions" | head -n 1)
    echo "\$conditions" | tail -n +2 | while read condition; do
        echo "\${condition}_vs_\${first_condition},\${condition},\${first_condition}" >> "${output_dir}/contrasts.csv"
    done
    
    # Download SRA data
    mkdir -p "${output_dir}/raw_reads"
    cd "${output_dir}/raw_reads"
    
    # Download SRA data
    while IFS=\$'\\t' read -r sra_id title info; do
        echo "Downloading \$sra_id..."
        # First prefetch the SRA data
        prefetch "\$sra_id"
        
        # Then use fasterq-dump with common parameters
        fasterq-dump \\
            --split-files \\
            --threads 4 \\
            --outdir . \\
            --skip-technical \\
            --read-filter pass \\
            --min-read-len 50 \\
            --max-read-len 1000 \\
            --include-technical \\
            --force \\
            "\$sra_id"
    done < "${output_dir}/sample_info.txt"
    """
} 