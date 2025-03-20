#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BAM_PREPROCESSING } from '../modules/bam_preprocessing'
include { PIRANHA_PEAK_CALLING } from '../modules/piranha_peak_calling'


// Parameters
params.gtf = file('gencode.v38.primary_assembly.annotation.gtf')
params.bam_dir = null
params.rmsk = null
params.outdir = 'results'
params.piranha_params = ''
params.genome_fasta = null

// Print usage message if required parameters are missing
if (!params.gtf || !params.bam_dir || !params.rmsk || !params.genome_fasta) {
    log.error """
    Required parameters missing!
    
    Usage:
    nextflow run main.nf --gtf GTF_FILE --bam_dir BAM_DIR --rmsk RMSK_BED --genome_fasta FASTA [options]
    
    Required arguments:
      --gtf          GTF file (e.g., gencode.v43.annotation.gtf)
      --bam_dir      Directory containing BAM files
      --rmsk         RepeatMasker BED file from UCSC
      --genome_fasta Reference genome FASTA file
    
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
    ${baseDir}/../bin/gtf_features_bedtools.sh -i $gtf -o .
    """
}



// Process to annotate peaks with genomic features
process ANNOTATE_FEATURES {
    publishDir "${params.outdir}/peak_annotations/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    path feature_dir
    
    output:
    tuple val(sample_id), path("*_annotated_peaks.bed"), path("*_feature_overlap_summary.csv"), path("*_feature_overlap_plot.pdf")
    
    script:
    """
    # Make script executable and copy it to current directory
    cp ${projectDir}/../bin/annotate_peaks.sh .
    chmod +x annotate_peaks.sh
    
    # Run script from current directory
    ./annotate_peaks.sh \\
        -p $peaks \\
        -f $feature_dir \\
        -o .
    """
}

// Process to annotate peaks with repeat elements
process ANNOTATE_REPEATS {
    publishDir "${params.outdir}/repeat_annotations/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    path rmsk
    
    output:
    tuple val(sample_id), path("*_rmsk_counts.bed"), path("*_rmsk_summary.csv"), path("*_rmsk_distribution.pdf")
    
    script:
    """
    # Make script executable and copy it to current directory
    cp ${projectDir}/bin/annotate_peaks_rmsk.sh .
    chmod +x annotate_peaks_rmsk.sh
    
    # Run script from current directory
    ./annotate_peaks_rmsk.sh \\
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
    path "summary_plots.pdf"
    
    script:
    """
    #!/usr/bin/env python3
    
    import glob
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from datetime import datetime
    
    # Set style for plots
    plt.style.use('seaborn')
    
    def create_feature_plot(df, title):
        plt.figure(figsize=(10, 6))
        sns.barplot(data=df, x='Feature', y='Count')
        plt.xticks(rotation=45, ha='right')
        plt.title(title)
        plt.tight_layout()
        return plt.gcf()
    
    def create_repeat_plot(df, title):
        plt.figure(figsize=(12, 6))
        sns.barplot(data=df, x='Repeat_Class', y='Count')
        plt.xticks(rotation=45, ha='right')
        plt.title(title)
        plt.tight_layout()
        return plt.gcf()
    
    # Initialize PDF for plots
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('summary_plots.pdf')
    
    # Process each sample
    sample_summaries = []
    for feature_summary in glob.glob('*/*_feature_overlap_summary.csv'):
        sample = feature_summary.split('/')[0]
        
        # Read feature statistics
        df_features = pd.read_csv(feature_summary)
        
        # Create feature distribution plot
        fig = create_feature_plot(df_features, f'Feature Distribution - {sample}')
        pdf.savefig(fig)
        plt.close()
        
        # Process repeat statistics if available
        rmsk_summary = glob.glob(f'*/{sample}/*_rmsk_summary.csv')
        if rmsk_summary:
            df_rmsk = pd.read_csv(rmsk_summary[0])
            fig = create_repeat_plot(df_rmsk, f'Repeat Element Distribution - {sample}')
            pdf.savefig(fig)
            plt.close()
        
        # Collect summary statistics
        summary = {
            'Sample': sample,
            'Total_Peaks': df_features['Count'].sum(),
            'Feature_Distribution': df_features.to_dict('records')
        }
        if rmsk_summary:
            summary['Repeat_Distribution'] = df_rmsk.to_dict('records')
        
        sample_summaries.append(summary)
    
    pdf.close()
    
    # Create HTML report
    html_content = []
    html_content.append('''
    <!DOCTYPE html>
    <html>
    <head>
        <title>Peak Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
            h1, h2 { color: #2c3e50; }
            .sample { margin: 20px 0; padding: 20px; background: #f7f9fc; border-radius: 5px; }
            table { border-collapse: collapse; width: 100%; margin: 10px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            .summary { background: #e8f4f8; padding: 15px; border-radius: 5px; margin: 10px 0; }
            .stats { display: flex; justify-content: space-around; flex-wrap: wrap; }
            .stat-box { 
                background: white; 
                padding: 15px; 
                margin: 10px; 
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                min-width: 200px;
            }
        </style>
    </head>
    <body>
        <h1>Peak Analysis Report</h1>
        <p><strong>Date:</strong> {}</p>
    '''.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
    # Add overall statistics
    total_peaks = sum(s['Total_Peaks'] for s in sample_summaries)
    html_content.append(f'''
        <div class="summary">
            <h2>Overall Statistics</h2>
            <div class="stats">
                <div class="stat-box">
                    <h3>Total Samples</h3>
                    <p>{len(sample_summaries)}</p>
                </div>
                <div class="stat-box">
                    <h3>Total Peaks</h3>
                    <p>{total_peaks:,}</p>
                </div>
                <div class="stat-box">
                    <h3>Average Peaks per Sample</h3>
                    <p>{total_peaks/len(sample_summaries):,.0f}</p>
                </div>
            </div>
        </div>
    ''')
    
    # Add sample-specific information
    for summary in sample_summaries:
        html_content.append(f'<div class="sample"><h2>Sample: {summary["Sample"]}</h2>')
        
        # Add feature statistics
        html_content.append('<h3>Feature Distribution</h3>')
        df_features = pd.DataFrame(summary['Feature_Distribution'])
        html_content.append(df_features.to_html(index=False))
        
        # Add repeat statistics if available
        if 'Repeat_Distribution' in summary:
            html_content.append('<h3>Repeat Element Distribution</h3>')
            df_rmsk = pd.DataFrame(summary['Repeat_Distribution'])
            html_content.append(df_rmsk.to_html(index=False))
        
        html_content.append('</div>')
    
    html_content.append('</body></html>')
    
    # Write HTML report
    with open('summary_report.html', 'w') as f:
        f.write('\\n'.join(html_content))
    
    # Create text report
    with open('summary_report.txt', 'w') as f:
        f.write('Peak Analysis Summary\\n')
        f.write('===================\\n')
        f.write(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\\n\\n')
        
        f.write(f'Total Samples: {len(sample_summaries)}\\n')
        f.write(f'Total Peaks: {total_peaks:,}\\n')
        f.write(f'Average Peaks per Sample: {total_peaks/len(sample_summaries):,.0f}\\n\\n')
        
        for summary in sample_summaries:
            f.write(f'Sample: {summary["Sample"]}\\n')
            f.write('Feature Distribution:\\n')
            df_features = pd.DataFrame(summary['Feature_Distribution'])
            f.write(df_features.to_string())
            f.write('\\n\\n')
            
            if 'Repeat_Distribution' in summary:
                f.write('Repeat Element Distribution:\\n')
                df_rmsk = pd.DataFrame(summary['Repeat_Distribution'])
                f.write(df_rmsk.to_string())
                f.write('\\n\\n')
    """
}

// Main workflow
workflow {
    // Input channel for BAM files
    bam_ch = Channel.fromPath("${params.bam_dir}/*.bam")
    // genome_fasta = Channel.fromPath(params.genome_fasta)

    // Run the workflow
    BAM_PREPROCESSING(bam_ch)
    PIRANHA_PEAK_CALLING(BAM_PREPROCESSING.out.processed_bam)
    
    // Extract features from GTF
    EXTRACT_FEATURES(params.gtf)
    
    // Annotate peaks with features
    ANNOTATE_FEATURES(
        PIRANHA_PEAK_CALLING.out.peaks,
        EXTRACT_FEATURES.out.feature_beds.collect()
    )
    
    // Annotate peaks with repeats
    ANNOTATE_REPEATS(
        PIRANHA_PEAK_CALLING.out.peaks,
        params.rmsk
    )
    
    // Generate final report
    GENERATE_REPORT(
        ANNOTATE_FEATURES.out.collect()
            .mix(ANNOTATE_REPEATS.out.collect())
            .collect()
    )
}