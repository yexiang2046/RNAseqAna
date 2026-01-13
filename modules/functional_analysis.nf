/*
 * process FUNCTIONAL_ANALYSIS
 * Performs comprehensive functional enrichment analysis on DEGs
 * - Gene Ontology (BP, MF, CC)
 * - KEGG pathways
 * - Reactome pathways
 * - GSEA with HALLMARK gene sets
 */
process FUNCTIONAL_ANALYSIS {
    tag "${comparison_name}"
    container 'xiang2019/rnaseq_renv:v1.0.0'
    publishDir "${params.outdir}/functional_analysis", mode: 'copy'

    input:
    tuple val(comparison_name), path(deg_file)

    output:
    path "functional_analysis_${comparison_name}/*", optional: true

    script:
    def species = params.species ?: 'human'
    """
    mkdir -p functional_analysis_${comparison_name}

    Rscript ${projectDir}/bin/functional_analysis.r \\
        ${deg_file} \\
        functional_analysis_${comparison_name} \\
        ${species} \\
        ${comparison_name}
    """
}
