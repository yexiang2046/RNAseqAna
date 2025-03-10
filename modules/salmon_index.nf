process SALMON_INDEX {
    publishDir "${params.outdir}/salmon_index", mode: 'copy'

    input:
    path transcriptome
    path genome

    output:
    path "salmon_index", emit: index

    script:
    """
    # Generate decoys file
    grep "^>" $genome | cut -d " " -f 1 | sed 's/>//g' > decoys.txt

    # Concatenate transcriptome and genome
    cat $transcriptome $genome > gentrome.fa

    # Create salmon index with decoys
    salmon index \
        -t gentrome.fa \
        -i salmon_index \
        -d decoys.txt \
        -k 31 \
        --threads $task.cpus
    """
} 