
import re
import subprocess
import os
from glob import glob
import itertools
#from bioinfodit.analysis import gff
import py_lib.utils as pu



subprocess.run(["bash", "source", "/home/xiang/miniconda3/bin/activate"])


# take fastq files and generate a dictionary of sample ID with respective fastq files
def get_samples(mode = "PE", fastq_file = fastq_files):
    while True:
        try:
            mode in ["PE", "SE"]
            break
        except ValueError:
            print("mode need to be PE or SE")
    samples = {}
    if mode == "PE":
        while isinstance(len(fastq_file)/2, int):
            try:
                for i in range(len(fastq_file)/2):
                    ID = os.path.basename(fastq_file[2*i]).split("_", 1)[0]
                    samples[ID] = [fastq_file[2*i], fastq_file[2*i+1]]
                break
            except ValueError:
                print("Make sure the fastq files are in pairs")
    else:
        for f in fastq_file:
            ID = os.path.basename(f).split("_", 1)[0]
            samples[ID] = f

    return samples



def trim_reads(files):
    



thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, Seq_mode, trimmed_output, fastq_dir, genome, annotation, counts_oufile = pu.get_para(configfile="config_template.txt")


fastq_files = get_fastq(configfile="config_template.txt")

samples = get_samples(mode = Seq_mode, fastq_file = fastq_files)


alignment_path = os.path.join(out_dir, alignment_dir)
trimmed_path = os.path.join(out_dir, trimmed_output)

if not os.path.exists(trimmed_path):
    os.mkdir(trimmed_path)

if not os.path.exists(alignment_path):
    os.mkdir(alignment_path)




if (Seq_mode == "PE"):
    for sample in samples.keys():
        subprocess.run(["./bash/step1_reads_trimming.sh", thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, Seq_mode, trimmed_path, samples[sample][0], samples[sample][1]])
else:
    for sample, fastq in samples:
        subprocess.run(["./bash/step1_reads_trimming.sh", thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, Seq_mode, trimmed_path, os.path.join(fastq_dir, fastq)])

for i in range(sample_num):
    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
#    print(file)
    trimmed_file1 = "".join([file, "_S1_L005_R1_001_val_1.fq.gz"])
    trimmed_files1 = os.path.join(trimmed_path, trimmed_file1)
    trimmed_file2 = "".join([file, "_S1_L005_R2_001_val_2.fq.gz"])
    trimmed_files2 = os.path.join(trimmed_path, trimmed_file2)
#    print(trimmed_files1, trimmed_files2)
    subprocess.run(["./bash/align.sh", thread_num, align_idx, out_dir, trimmed_files1, trimmed_files2, alignment_path])


for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    bamfile = "".join([alignment_path, "/", file, "_S1_L005_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam"])
    outputfile = "".join([alignment_path, "/", file, ".filtered.bam"])
    index_bam(bamfile, thread = thread_num)
    filter_bam(bamfile = bamfile, thread = thread_num, outdir = alignment_path, output_file = outputfile, genome = genome)

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
