#!/usr/bin/env python3
import re
from sys import stderr, exit
from argparse import ArgumentParser, FileType

def extract_genes_gtf(input_gff_path, gene_list_path, output_gff_path):
    print("Input GTF File:", input_gff_path)
    print("Gene List File:", gene_list_path)
    print("Output GTF File:", output_gff_path)
    with open(output_gff_path, "w") as outfile:
        gene_names = set()
        with open(gene_list_path, 'r') as gene_file:
            for gene in gene_file:
                gene_names.add(gene.strip())

        with open(input_gff_path, 'r') as features:
            for feature in features:
                for gene in gene_names:
                    regex = 'gene_name "' + gene + '";'
                    if re.search(regex, feature):
                        outfile.write(feature)
                        break  # Stop searching for this feature in the gene list


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract GTF item from GTF file by gene names')
    parser.add_argument('input_gtf_file',
            type=str,
            help='input GTF file')
    parser.add_argument('gene_list',
            type=str,
            help='input gene list file, one gene per line')
    parser.add_argument('output_gtf_file',
            type=str,
            help='output GTF file')
    parser.add_argument('-v', '--verbose',
            dest='verbose',
            action='store_true',
            help='also print some statistics to stderr')

args = parser.parse_args()
extract_genes_gtf(args.input_gtf_file, args.gene_list, args.output_gtf_file)

