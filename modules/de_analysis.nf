/*
 * define the DE_ANALYSIS process that performs differential expression analysis
 * given the counts file and the metadata file
 */
process DE_ANALYSIS {
    debug true
    container 'xiang2019/rnaseq_renv:v1.0.1'
    publishDir "${params.outdir}", mode:'copy'
    input:
    path counts_file
    path metadata_file
    path gtffile
    val sp
    
    output:
    path "de_results/*.csv"
    path "de_results/PCA_plot.png"

    script:
    """
    mkdir -p de_results
    edger.r -c $counts_file -g $gtffile  -m $metadata_file -o "de_results" -s $sp
    """
}
