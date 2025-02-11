/*
 * define the DE_ANALYSIS process that performs differential expression analysis
 * given the counts file and the metadata file
 */
process DE_ANALYSIS {
    debug true
    container 'xiang2019/rnaseq_renv:v1.0.0'
	publishDir "${projectDir}/de_results", mode:'copy'
    input:
    path counts_file
    path metadata_file

    output:
    path "*.csv"
    path "PCA_plot.png"
    path "DEG_barplot_*.png"

    script:
    """
    Rscript path/to/differential_expression_analysis.R --counts $counts_file --metadata $metadata_file --output "de_results"
    """
}