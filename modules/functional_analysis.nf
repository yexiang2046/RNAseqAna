/*
 * process FUNCTIONAL_ANALYSIS
 * performs functional analysis on the differentially expressed genes
 * uses the clusterProfiler and fgsea packages
 * uses the msigdb database
 */
process FUNCTIONAL_ANALYSIS {
    debug true
    container 'xiang2019/rnaseq_renv:v1.0.0'
    publishDir "${projectDir}/functional_analysis_results", mode:'copy'
    input:
    path de_results_files from de_results_ch

    output:
    path "functional_analysis_results/*"

    script:
    """
    mkdir -p functional_analysis_results
    Rscript bin/funtional_analysis.r
    """
}