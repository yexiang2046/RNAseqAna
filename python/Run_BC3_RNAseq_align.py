
import re
import subprocess
import os
from glob import glob
import itertools
#from bioinfodit.analysis import gff
import py_lib.utils as pu



# subprocess.run(["bash", "source", "/home/xiang/miniconda3/bin/activate"])




    



thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, Seq_mode, trimmed_output, fastq_dir, genome, annotation, counts_oufile = pu.get_para(configfile="config_BC3_RNAseq.txt")


fastq_files = pu.get_fastq(configfile="config_BC3_RNAseq.txt")

print(fastq_files)

alignment_path = os.path.join(out_dir, alignment_dir)
trimmed_path = os.path.join(out_dir, trimmed_output)

if not os.path.exists(trimmed_path):
    os.mkdir(trimmed_path)

if not os.path.exists(alignment_path):
    os.mkdir(alignment_path)

sample_num = int(len(fastq_files)/2)

for i in range(sample_num):
    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
#    print(file)
    file1 = "".join([file, "_S1_L005_R1_001.fastq.gz"])
    files1 = os.path.join(fastq_dir, file1)
    file2 = "".join([file, "_S1_L005_R2_001.fastq.gz"])
    files2 = os.path.join(fastq_dir, file2)
#    print(trimmed_files1, trimmed_files2)
    subprocess.run(["./bash/align.sh", thread_num, align_idx, out_dir, files1, files2, alignment_path])


for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    bamfile = "".join([alignment_path, "/", file, "_S1_L005_R1_001.fastq.gzAligned.sortedByCoord.out.bam"])
    outputfile = "".join([alignment_path, "/", file, ".filtered.bam"])
    pu.index_bam(bamfile, thread = thread_num)
    pu.filter_bam(bamfile = bamfile, thread = thread_num, outdir = alignment_path, output_file = outputfile, genome = genome)

for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
    index_bam(filtered_bam, thread = thread_num)

bw_dir = os.path.join(out_dir, "bw_files")

if not os.path.exists(bw_dir):
    os.mkdir(bw_dir)


for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
    bigwig_file = "".join([bw_dir, "/", file, ".bw"])

    bam_to_bw(bamfile = filtered_bam, thread=thread_num, output_file=bigwig_file)

bam_list =[]
for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
    bam_list.append(filtered_bam)


GTF = annotation
GTF_file = os.path.join(out_dir, GTF)
count_file = counts_outfile
count_path_file = os.path.join(out_dir, count_file)
#print(count_file)
#print(count_path_file)
#print(bam_list)
count_features(bamfiles = bam_list, gtf=GTF_file, thread=thread_num, count_outfile=count_path_file)


## Generate report files
subprocess.run(["multiqc", out_dir, "-o", out_dir])
