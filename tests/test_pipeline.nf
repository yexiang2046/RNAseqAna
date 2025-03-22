#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import main workflow
include { run_piranha } from '../workflows/run_piranha'

// Test parameters
params.test_data_dir = "$baseDir/data"
params.outdir = "test_results"

// Test data paths
def test_bam = file("${params.test_data_dir}/bamfiles/test.bam")
def test_gtf = file("${params.test_data_dir}/gencode.v38.primary_assembly.annotation.gtf")
def test_rmsk = file("${params.test_data_dir}/all_rmsk_hg38.bed")
def test_genome = file("${params.test_data_dir}/GRCh38.primary_assembly.genome.fa")

// Test process to create mock BAM file
process CREATE_TEST_BAM {
    output:
    path "test.bam"
    
    script:
    """
    # Create a small BAM file with simulated header and read
    echo -e "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:248956422
@PG\tID:test\tPN:test
test_read1\t0\tchr1\t100\t255\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" > test.sam
    samtools view -b test.sam > test.bam
    samtools index test.bam
    """
}

// Test process to create mock GTF file
process CREATE_TEST_GTF {
    output:
    path "test.gtf"
    
    script:
    """
    # Create a small GTF file with test data
    cat > test.gtf << 'EOF'
    chr1    HAVANA  gene    100     200     .       +       .       gene_id "TEST1"; gene_name "TestGene1"; gene_type "protein_coding";
    chr1    HAVANA  exon    100     150     .       +       .       gene_id "TEST1"; transcript_id "TEST1.1"; exon_number "1";
    chr1    HAVANA  exon    160     200     .       +       .       gene_id "TEST1"; transcript_id "TEST1.1"; exon_number "2";
    chr1    HAVANA  CDS     120     180     .       +       0       gene_id "TEST1"; transcript_id "TEST1.1"; protein_id "TEST1.1";
    EOF
    """
}

// Test process to create mock RMSK file
process CREATE_TEST_RMSK {
    output:
    path "test_rmsk.bed"
    
    script:
    """
    # Create a small RMSK file with test data
    cat > test_rmsk.bed << 'EOF'
    chr1    100     150     L1      LINE    +       1000    0.95
    chr1    200     250     Alu     SINE    -       800     0.98
    EOF
    """
}

// Test process to validate BAM preprocessing
process VALIDATE_BAM_PREPROCESSING {
    input:
    path bam
    
    output:
    path "validation_report.txt"
    
    script:
    """
    # Validate BAM file
    echo "Validating BAM file: $bam" > validation_report.txt
    samtools view -H $bam >> validation_report.txt
    echo "Number of reads: \$(samtools view -c $bam)" >> validation_report.txt
    echo "BAM index exists: \$(test -f ${bam}.bai && echo 'Yes' || echo 'No')" >> validation_report.txt
    """
}

// Test process to validate peak calling
process VALIDATE_PEAK_CALLING {
    input:
    path peaks
    
    output:
    path "validation_report.txt"
    
    script:
    """
    # Validate peak file
    echo "Validating peak file: $peaks" > validation_report.txt
    echo "Number of peaks: \$(wc -l < $peaks)" >> validation_report.txt
    echo "File format check:" >> validation_report.txt
    head -n 1 $peaks >> validation_report.txt
    """
}

// Test process to validate feature annotation
process VALIDATE_FEATURE_ANNOTATION {
    input:
    path annotated_peaks
    path summary_csv
    path plot_pdf
    
    output:
    path "validation_report.txt"
    
    script:
    """
    # Validate annotation outputs
    echo "Validating feature annotation outputs" > validation_report.txt
    echo "Annotated peaks file exists: \$(test -f $annotated_peaks && echo 'Yes' || echo 'No')" >> validation_report.txt
    echo "Summary CSV exists: \$(test -f $summary_csv && echo 'Yes' || echo 'No')" >> validation_report.txt
    echo "Plot PDF exists: \$(test -f $plot_pdf && echo 'Yes' || echo 'No')" >> validation_report.txt
    """
}

// Test workflow
workflow {
    // Create test data
    template_bam = Channel.fromPath("${params.test_data_dir}/bamfiles/test.bam")
    test_bam = CREATE_TEST_BAM()
    test_gtf = CREATE_TEST_GTF()
    test_rmsk = CREATE_TEST_RMSK()
    
    // Run main workflow with test data
    main_workflow = run_piranha(
        test_bam,
        test_gtf,
        test_rmsk,
        Channel.fromPath("${params.test_data_dir}/GRCh38.primary_assembly.genome.fa"),
        params.outdir,
        ''
    )
    
    // Validate outputs
    VALIDATE_BAM_PREPROCESSING(main_workflow.processed_bam)
    VALIDATE_PEAK_CALLING(main_workflow.peaks_bed)
    VALIDATE_FEATURE_ANNOTATION(
        main_workflow.annotated_peaks,
        main_workflow.feature_summary,
        main_workflow.feature_plot
    )
    
    // Print test results
    main_workflow.report_html.view { "Generated HTML report: $it" }
    main_workflow.report_txt.view { "Generated text report: $it" }
    main_workflow.report_plots.view { "Generated plots: $it" }
} 