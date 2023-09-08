#!/bin/bash

#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0001.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0001.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0002.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0002.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0003.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0003.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0004.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0004.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0005.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0005.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0006.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0006.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0007.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0007.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0008.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0008.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0009.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0009.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0010.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0010.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0011.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0011.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-YZ-0012.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0012.unique.bed

#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0001.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0001.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0002.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0002.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0003.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0003.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0004.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0004.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0005.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0005.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0006.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0006.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0007.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0007.unique.bed
#bedtools bamtobed -i /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-YZ-0008.unique.op.bam > /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0008.unique.bed


# using docker image nfcore/clipseq
Piranha -o wt_inf6h_piranha_outa /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0003.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0001.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0002.unique.bed -z 200 -l -s
Piranha -o wt_inf6h_piranha_outb /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0004.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0001.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0002.unique.bed -z 200 -l -s

Piranha -o ko_inf6h_piranha_outa /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0007.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0003.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0004.unique.bed -z 200 -l -s
Piranha -o ko_inf6h_piranha_outb /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0008.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0003.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9814-yz-0004.unique.bed -z 200 -l -s



Piranha -o wt_uninf_piranha_outa /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0007.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0001.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0002.unique.bed -z 200 -l -s
Piranha -o wt_uninf_piranha_outb /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0008.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0001.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0002.unique.bed -z 200 -l -s

Piranha -o wt_inf24h_piranha_outa /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0009.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0003.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0004.unique.bed -z 200 -l -s
Piranha -o wt_inf24h_piranha_outb /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0010.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0003.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0004.unique.bed -z 200 -l -s

Piranha -o ko_inf24h_piranha_outa /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0011.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0005.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0006.unique.bed -z 200 -l -s
Piranha -o ko_inf24h_piranha_outb /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0012.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0005.unique.bed /home/xiang/New_disk2/vacv_j2_ip_yz/aligned/9221-yz-0006.unique.bed -z 200 -l -s


# annotate Piranha peaks by bedtools intersect
# bedtools intersect -wo -a WT_24h_shared.bed -b hg38_rmsk_gencodev43.gtf > WT_24h_shared_annotation.bed
# bedtools intersect -wo -a KO_24h_shared.bed -b hg38_rmsk_gencodev43.gtf > KO_24h_shared_annotation.bed
# bedtools intersect -wo -a KO_6h_shared.bed -b hg38_rmsk_gencodev43.gtf > KO_6h_shared_annotation.bed
# bedtools intersect -wo -a WT_6h_shared.bed -b hg38_rmsk_gencodev43.gtf > WT_6h_shared_annotation.bed
# bedtools intersect -wo -a WT_uninf_shared.bed -b hg38_rmsk_gencodev43.gtf > WT_uninf_shared_annotation.bed
# bedtools intersect -wo -a WT_uninf_shared.bed -b hg38_rmsk_gencodev43.gtf > WT_uninf_shared_annotation.bed
# bedtools intersect -wo -a WT_6h_shared.bed -b hg38_rmsk_gencodev43.gtf > WT_6h_shared_annotation.bed
# bedtools intersect -wo -a KO_6h_shared.bed -b hg38_rmsk_gencodev43.gtf > KO_6h_shared_annotation.bed
# bedtools intersect -wo -a WT_24h_shared.bed -b hg38_rmsk_gencodev43.gtf > WT_24h_shared_annotation.bed
# bedtools intersect -wo -a KO_24h_shared.bed -b hg38_rmsk_gencodev43.gtf > KO_24h_shared_annotation.bed