#!/bin/bash

source ~/miniconda3/bin/activate
conda activate fastp

# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-1_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-1.trimmed.CGATGT.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-1.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-2_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-2.trimmed.TGACCA.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-2.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-3_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-3.trimmed.CACTGT.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-3.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-4_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-4.trimmed.ATTGGC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-4.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-5_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-5.trimmed.CAGATC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-5.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-6_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-6.trimmed.TACAAG.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-6.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-7_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-7.trimmed.CCGTCC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-7.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-8_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-8.trimmed.TTTCAC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-8.fastq.html
# fastp -i ~/New_disk1/PROseq_PEL/OG/8199-XY-9_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-9.trimmed.CGATGT.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-9.fastq.html
 


fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-1_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-1.trimmed.CGATGT.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-1.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-2_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-2.trimmed.TGACCA.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-2.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-3_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-3.trimmed.CACTGT.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-3.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-4_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-4.trimmed.ATTGGC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-4.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-5_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-5.trimmed.CAGATC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-5.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-6_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-6.trimmed.TACAAG.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-6.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-7_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-7.trimmed.CCGTCC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-7.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-8_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-8.trimmed.TTTCAC.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-8.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/8199-XY-9_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/8199-XY-9.trimmed.CGATGT.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/8199-XY-9.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
 
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-1_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-1.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-1.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-2_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-2.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-2.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-3_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-3.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-3.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-4_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-4.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-4.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-5_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-5.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-5.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-6_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-6.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-6.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-7_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-7.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-7.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-8_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-8.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-8.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-9_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-9.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-9.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa
fastp -p -y -Y 50 -i ~/New_disk1/PROseq_PEL/OG/7667-XY-10_S1_L005_R1_001.fastq.gz -o ~/New_disk1/PROseq_PEL/OG/7667-XY-10.trimmed.fastq.gz -h ~/New_disk1/PROseq_PEL/OG/7667-XY-10.fastq.html -q 30 --adapter_fasta PROseq_adaptor.fa