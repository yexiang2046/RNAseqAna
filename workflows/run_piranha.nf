#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PIRANHA_PEAK_CALLING } from '../modules/piranha_peak_calling.nf'


// Parameters
params.gtf = 'gencode.v38.primary_assembly.annotation.gtf'
params.bam_dir = null
params.rmsk = null
params.outdir = 'results'
params.piranha_params = ''

// Print usage message if required parameters are missing
if (!params.gtf || !params.bam_dir || !params.rmsk) {
    log.error """
    Required parameters missing!
    
    Usage:
    nextflow run main.nf --gtf GTF_FILE --bam_dir BAM_DIR --rmsk RMSK_BED [options]
    
    Required arguments:
      --gtf         GTF file (e.g., gencode.v43.annotation.gtf)
      --bam_dir     Directory containing BAM files
      --rmsk        RepeatMasker BED file from UCSC
    
    Optional arguments:
      --outdir            Output directory (default: 'results')
      --piranha_params    Additional parameters for Piranha
    """
    exit 1
}

// Process to extract genomic features from GTF
process EXTRACT_FEATURES {
    publishDir "${params.outdir}/features_bed", mode: 'copy'
    
    input:
    path gtf
    
    output:
    path "*.bed", emit: feature_beds
    
    script:
    """
    gtf_features_bedtools.sh -i $gtf -o .
    """
}



// Process to annotate peaks with genomic features
process ANNOTATE_FEATURES {
    publishDir "${params.outdir}/peak_annotations/\${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    path feature_dir
    
    output:
    tuple val(sample_id), path("*_annotated_peaks.bed"), path("*_feature_overlap_summary.csv"), path("*_feature_overlap_plot.pdf")
    
    script:
    """
    annotate_peaks.sh \\
        -p $peaks \\
        -f $feature_dir \\
        -o .
    """
}

// Process to annotate peaks with repeat elements
process ANNOTATE_REPEATS {
    publishDir "${params.outdir}/repeat_annotations/\${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    path rmsk
    
    output:
    tuple val(sample_id), path("*_rmsk_counts.bed"), path("*_rmsk_summary.csv"), path("*_rmsk_distribution.pdf")
    
    script:
    """
    annotate_peaks_rmsk.sh \\
        -p $peaks \\
        -r $rmsk \\
        -o .
    """
}

// Process to generate final report
process GENERATE_REPORT {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path '*'
    
    output:
    path "summary_report.html"
    path "summary_report.txt"
    
    script:
    """
    #!/usr/bin/env python3
    import glob
    import pandas as pd
    from datetime import datetime
    
    # Create HTML report
    html_content = []
    html_content.append('''
    <!DOCTYPE html>
    <html>
    <head>
        <title>Piranha Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            h1, h2 { color: #2c3e50; }
            .sample { margin: 20px 0; padding: 20px; background: #f7f9fc; border-radius: 5px; }
            table { border-collapse: collapse; width: 100%; margin: 10px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
        </style>
    </head>
    <body>
        <h1>Piranha Analysis Report</h1>
        <p><strong>Date:</strong> {}</p>
    '''.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
    # Process each sample
    for feature_summary in glob.glob('*/*_feature_overlap_summary.csv'):
        sample = feature_summary.split('/')[0]
        html_content.append(f'<div class="sample"><h2>Sample: {sample}</h2>')
        
        # Add feature statistics
        df = pd.read_csv(feature_summary)
        html_content.append('<h3>Feature Overlap Statistics</h3>')
        html_content.append(df.to_html(index=False))
        
        # Add repeat statistics if available
        rmsk_summary = glob.glob(f'*/{sample}/*_rmsk_summary.csv')
        if rmsk_summary:
            df_rmsk = pd.read_csv(rmsk_summary[0])
            html_content.append('<h3>Repeat Element Statistics</h3>')
            html_content.append(df_rmsk.to_html(index=False))
        
        html_content.append('</div>')
    
    html_content.append('</body></html>')
    
    # Write HTML report
    with open('summary_report.html', 'w') as f:
        f.write('\\n'.join(html_content))
    
    # Create text report
    with open('summary_report.txt', 'w') as f:
        f.write('Piranha Analysis Summary\\n')
        f.write('=======================\\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\\n\\n')
        
        for feature_summary in glob.glob('*/*_feature_overlap_summary.csv'):
            sample = feature_summary.split('/')[0]
            f.write(f'Sample: {sample}\\n')
            f.write('Feature Overlap Statistics:\\n')
            df = pd.read_csv(feature_summary)
            f.write(df.to_string())
            f.write('\\n\\n')
            
            rmsk_summary = glob.glob(f'*/{sample}/*_rmsk_summary.csv')
            if rmsk_summary:
                f.write('Repeat Element Statistics:\\n')
                df_rmsk = pd.read_csv(rmsk_summary[0])
                f.write(df_rmsk.to_string())
                f.write('\\n\\n')
    """
}

// Main workflow
workflow {
    // Get input files
    gtf_file = file(params.gtf)
    rmsk_file = file(params.rmsk)
    bam_files = Channel.fromPath("${params.bam_dir}/*.bam")
    
    // Extract features from GTF
    EXTRACT_FEATURES(gtf_file)
    
    // Run Piranha on each BAM file
    PIRANHA_PEAK_CALLING(bam_files)
    
    // Annotate peaks with features
    ANNOTATE_FEATURES(
        PIRANHA_PEAK_CALLING.out.peaks,
        EXTRACT_FEATURES.out.feature_beds.collect()
    )
    
    // Annotate peaks with repeats
    ANNOTATE_REPEATS(
        PIRANHA_PEAK_CALLING.out.peaks,
        rmsk_file
    )
    
    // Generate final report
    GENERATE_REPORT(
        ANNOTATE_FEATURES.out.collect()
            .mix(ANNOTATE_REPEATS.out.collect())
            .collect()
    )
}