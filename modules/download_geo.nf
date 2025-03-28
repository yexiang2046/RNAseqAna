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
    
    # Download GEO metadata
    GEO_URL="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${geo_accession}&targ=self&form=text&view=full"
    curl -s "\$GEO_URL" > "${output_dir}/${geo_accession}_metadata.txt"
    
    # Extract SRA IDs and sample information
    awk '
        /^!Sample_geo_accession/ { split(\$0, a, "="); gsub(/^[ \\t]+/, "", a[2]); sra[NR] = a[2] }
        /^!Sample_characteristics_ch1/ { split(\$0, a, "="); gsub(/^[ \\t]+/, "", a[2]); info[NR] = a[2] }
        /^!Sample_title/ { split(\$0, a, "="); gsub(/^[ \\t]+/, "", a[2]); title[NR] = a[2] }
        END {
            for (i in sra) {
                if (sra[i] != "") {
                    print sra[i] "\\t" title[i] "\\t" info[i]
                }
            }
        }
    ' "${output_dir}/${geo_accession}_metadata.txt" > "${output_dir}/sample_info.txt"
    
    # Create metadata.csv for edgeR
    echo "sample,condition" > "${output_dir}/metadata.csv"
    awk -F'\\t' '
        {
            # Extract condition from characteristics (assuming format "condition: value")
            split(\$3, a, ":")
            if (length(a) > 1) {
                gsub(/^[ \\t]+/, "", a[2])
                print \$1 "," a[2]
            } else {
                # If no characteristics, use sample title
                print \$1 "," \$2
            }
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