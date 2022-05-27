
import re
import subprocess
import os
from glob import glob
import itertools
import py_lib.utils as pu
#from bioinfodit.analysis import gff





Para = pu.get_para(configfile="config_PCIF_exoseq.txt")
thread_num = Para[0]
align_idx = Para[1]
out_dir = Para[2]
trimmer_dir = Para[3]
alignment_dir = Para[4]
pair_end = Para[5]
trimmed_output = Para[6]
fastq_dir = Para[7]
genome = Para[8]
annotation = Para[9]
counts_outfile = Para[10]


######################################
# change the config file accorddingly
fastq_files = pu.get_fastq(configfile="config_PCIF_exoseq.txt")

samples = []

for i in range(int(len(fastq_files)/2)):
    samples.append(os.path.basename(fastq_files[2*i]).split('_', 1)[0])


sample_num = len(samples)

alignment_path = os.path.join(out_dir, alignment_dir)


if not os.path.exists(alignment_path):
    os.mkdir(alignment_path)


if (pair_end == "Yes"):


    for i in range(sample_num):
        #print(file)
        read1_file = fastq_files[2*i]
        read1_file_Wpath = os.path.join(fastq_dir, read1_file)
        read2_file = fastq_files[2*i+1]
        read2_file_Wpath = os.path.join(fastq_dir, read2_file)
        output_file = "".join([samples[i], ".sam"])
        output_file_Wpath = os.path.join(alignment_path, output_file)
        print("processing ", samples[i])
        # pu.hisat2_align(thread_num = thread_num, idx = align_idx, reads_1 = read1_file_Wpath, reads_2 = read2_file_Wpath, output_sam = output_file_Wpath)

elif (pair_end == "No"):
    # not working
    print("single end alignment is not supported yet!")
else:
    print("Answer to pair_end should be Yes or No")

print(genome)
for i in range(sample_num):

        sam_file = "".join([samples[i], ".sam"])
        sam_file_Wpath = os.path.join(alignment_path, sam_file)
        bam_file_Wpath = os.path.join(alignment_path, "".join([samples[i], ".bam"]))
        bam_file_sort_Wpath = "".join([bam_file_Wpath, ".srt"])
        # pu.sam2bam_q10(sam = sam_file_Wpath, bam = bam_file_Wpath)
        # pu.sortbam(bam=bam_file_Wpath, sortedBam=bam_file_sort_Wpath, thread_num=thread_num)
        filtered_bam = "".join([alignment_path, "/", samples[i], ".filtered.bam"])
        # pu.index_bam(bam_file_sort_Wpath, thread = thread_num)
        pu.filter_bam(bamfile = bam_file_sort_Wpath, thread = thread_num, outdir = alignment_path, output_file = filtered_bam, genome = genome)
        pu.index_bam(filtered_bam, thread = thread_num)

bw_dir = os.path.join(out_dir, "bw_files")


if not os.path.exists(bw_dir):
    os.mkdir(bw_dir)


for i in range(sample_num):

    bamfile = os.path.join(alignment_path, "".join([samples[i], ".filtered.bam"]))
    bigwig_file = "".join([bw_dir, "/", samples[i], ".bw"])
    # pu.bam_to_bw(bamfile = bamfile, thread=thread_num, output_file=bigwig_file, norm_method = "BPM")

bam_list =[]
for i in range(sample_num):

        bamfile = "".join([alignment_path, "/", samples[i], ".filtered.bam"])
        bam_list.append(bamfile)


GTF = annotation
GTF_file = os.path.join(out_dir, GTF)
count_file = counts_outfile
count_path_file = os.path.join(out_dir, count_file)
#print(count_file)
#print(count_path_file)
#print(bam_list)
# pu.count_features(bamfiles = bam_list, gtf=GTF_file, thread=thread_num, count_outfile=count_path_file)


## Generate report files
# subprocess.run(["multiqc", out_dir, "-o", out_dir])
