#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import main workflow
include { run_piranha } from '../workflows/run_piranha'

// Test parameters
params.test_data_dir = "$baseDir/data"
params.outdir = "test_results"

// Test process to create mock BAM file
process CREATE_TEST_BAM {
    output:
    path "test.bam", emit: bam
    path "test.bam.bai", emit: bai
    
    script:
    """
    # Create SAM header
    cat << EOF > test.sam
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@PG	ID:test	PN:test
EOF
    
    # Define sequence and quality strings
    seq="ATCGATCGATATCGATCGAT"
    qual="HHHHHHHHHHHHHHHHHHHH"  # Higher quality scores
    
    # Create paired-end reads for peak regions
    # First peak region (100-120)
    for pos in {100..120..2}; do
        mate1_pos=\$pos
        mate2_pos=\$((pos + 40))  # 40bp insert size
        read_name="read_\${pos}"
        
        # Mate 1 (forward)
        echo -e "\${read_name}\t99\tchr1\t\${mate1_pos}\t60\t20M\t=\t\${mate2_pos}\t60\t\${seq}\t\${qual}" >> test.sam
        
        # Mate 2 (reverse)
        echo -e "\${read_name}\t147\tchr1\t\${mate2_pos}\t60\t20M\t=\t\${mate1_pos}\t-60\t\${seq}\t\${qual}" >> test.sam
    done
    
    # Second peak region (200-220)
    for pos in {200..220..2}; do
        mate1_pos=\$pos
        mate2_pos=\$((pos + 40))
        read_name="read_\${pos}"
        
        # Mate 1 (forward)
        echo -e "\${read_name}\t99\tchr1\t\${mate1_pos}\t60\t20M\t=\t\${mate2_pos}\t60\t\${seq}\t\${qual}" >> test.sam
        
        # Mate 2 (reverse)
        echo -e "\${read_name}\t147\tchr1\t\${mate2_pos}\t60\t20M\t=\t\${mate1_pos}\t-60\t\${seq}\t\${qual}" >> test.sam
    done
    
    # Convert to BAM and sort by coordinate
    samtools view -b test.sam | samtools sort -o test.bam
    samtools index test.bam
    
    # Validate BAM file
    echo "=== BAM validation ===" 
    samtools validate test.bam || true
    echo "=== BAM statistics ==="
    samtools flagstat test.bam || true
    echo "=== First few reads ==="
    samtools view test.bam | head -n 4 || true
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

// Test process to create mock FASTA file
process CREATE_TEST_FASTA {
    output:
    path "test.fa"
    
    script:
    """
    # Create a small FASTA file with test sequence
    cat > test.fa << 'EOF'
>chr1 GRCh38.p13 Primary Assembly
AGCTATCACAGATCGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
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
    CREATE_TEST_BAM()
    test_gtf = CREATE_TEST_GTF()
    test_rmsk = CREATE_TEST_RMSK()
    test_fa = CREATE_TEST_FASTA()
    
    // Create channels for the test data
    bam_ch = CREATE_TEST_BAM.out.bam
    gtf_ch = test_gtf
    rmsk_ch = test_rmsk
    genome_ch = test_fa
    
    // Run main workflow with test data
    main_workflow = run_piranha(
        bam_ch,
        gtf_ch,
        rmsk_ch,
        genome_ch,
        params.outdir,
        ''
    )
    
    // Validate outputs
    main_workflow.processed_bam
        .map { tuple -> 
            println "Debug - Full tuple: ${tuple}"
            // Extract just the BAM file path from the tuple
            def bam_path = tuple[1]  // Get the second element (index 1) which is the path
            println "Debug - Extracted BAM path: ${bam_path}"
            return bam_path
        }
        .set { bam_for_validation }
    
    VALIDATE_BAM_PREPROCESSING(bam_for_validation)
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