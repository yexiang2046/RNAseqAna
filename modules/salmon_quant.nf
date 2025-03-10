process SALMON_QUANT {
    tag "$sample_id"
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant_results
    path "${sample_id}/quant.sf", emit: quant_sf
    path "${sample_id}/lib_format_counts.json", emit: lib_format
    path "${sample_id}/aux_info/meta_info.json", emit: meta_info
    path "${sample_id}"

    script:
    def single_end = reads instanceof Path
    
    if (single_end) {
        """
        salmon quant \
            -i $index \
            -l A \
            -r $reads \
            --validateMappings \
            --gcBias \
            --seqBias \
            -o ${sample_id} \
            --threads $task.cpus
        """
    } else {
        """
        salmon quant \
            -i $index \
            -l A \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            --validateMappings \
            --gcBias \
            --seqBias \
            -o ${sample_id} \
            --threads $task.cpus
        """
    }
} 