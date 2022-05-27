
import re
import subprocess
import os
from glob import glob
import itertools
import py_lib.utils as pu
#from bioinfodit.analysis import gff





Para = pu.get_para(configfile="config_myeloid_RNAseq.txt")
thread_num = Para[0]
align_idx = Para[1]
out_dir = Para[2]
trimmer_dir = Para[3]
alignment_dir = Para[4]
pair_end = Para[5]
trimmed_output = Para[6]
fastq_dir = Para[7]
genome = Para[8]
# annotation = Para[9]
# counts_outfile = Para[10]


######################################
# change the config file accorddingly
fastq_files = pu.get_fastq(configfile="config_myeloid_RNAseq.txt")

samples = []

for f in fastq_files:
    samples.append(os.path.basename(f).split('_', 1)[0])

alignment_path = os.path.join(out_dir, alignment_dir)
trimmed_path = os.path.join(out_dir, trimmed_output)

if not os.path.exists(trimmed_path):
    os.mkdir(trimmed_path)

if not os.path.exists(alignment_path):
    os.mkdir(alignment_path)


print(Para)
print(fastq_files)
sample_num = []
if (pair_end == "Yes"):
    sample_num = int(len(fastq_files)/2)
    for i in range(sample_num):
        print(i)
        subprocess.run(["./bash/step1_reads_trimming.sh", thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, pair_end, trimmed_path, fastq_files[i*2], fastq_files[i*2+1]])
else:
    sample_num = len(fastq_files)
    for i in range(len(fastq_files)):
        subprocess.run(["./bash/step1_reads_trimming.sh", thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, pair_end, trimmed_path, os.path.join(fastq_dir, fastq_files[i])])
        print(os.path.join(fastq_dir, fastq_files[i]))

if (pair_end == "Yes"):


    for i in range(sample_num):
        file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
        #print(file)
        trimmed_file1 = "".join([file, "_S1_L005_R1_001_val_1.fq.gz"])
        trimmed_files1 = os.path.join(trimmed_path, trimmed_file1)
        trimmed_file2 = "".join([file, "_S1_L005_R2_001_val_2.fq.gz"])
        trimmed_files2 = os.path.join(trimmed_path, trimmed_file2)
        #print(trimmed_files1, trimmed_files2)
        subprocess.run(["./bash/align.sh", thread_num, align_idx, out_dir, trimmed_files1, trimmed_files2, alignment_path])
elif (pair_end == "No"):
    for i in range(sample_num):
        file = os.path.basename(fastq_files[i]).split('_', 1)[0]
        trimmed_file = "".join([file, "_S1_L005_R1_001_trimmed.fq.gz"])
        trimmed_files = os.path.join(trimmed_path, trimmed_file)
        subprocess.run(["./bash/align.sh", thread_num, align_idx, out_dir, trimmed_files, "No", alignment_path])

else:
    print("Answer to pair_end should be Yes or No")

for i in range(sample_num):

        file = os.path.basename(fastq_files[i]).split('_', 1)[0]
        filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
        bamfile = "".join([alignment_path, "/", file, "_S1_L005_R1_001_trimmed.fq.gzAligned.sortedByCoord.out.bam"])
        pu.index_bam(bamfile, thread = thread_num)
        pu.filter_bam(bamfile = bamfile, thread = thread_num, outdir = alignment_path, output_file = filtered_bam, genome = genome)
        pu.index_bam(filtered_bam, thread = thread_num)

bw_dir = os.path.join(out_dir, "bw_files")


if not os.path.exists(bw_dir):
    os.mkdir(bw_dir)


for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
    bigwig_file = "".join([bw_dir, "/", file, ".bw"])

    #bam_to_bw(bamfile = filtered_bam, thread=thread_num, output_file=bigwig_file)

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
