#!/usr/bin/env python3

from sys import stderr, exit
from argparse import ArgumentParser, FileType
import numpy as np
import pandas as pd

def extend_splice_sites(splice_site_file, out_bedfile):
    # postitions_file format with four column, chr, left, right, strand
    
    Five_prime_bed = out_bedfile + ".5p.bed"
    Three_prime_bed = out_bedfile + ".3p.bed"
   
    positions = pd.read_table(splice_site_file, header=None)
    positions = positions.rename(columns={0:"Chr", 1:"Five", 2:"Three", 3:"Strand"})
    
    for ind in positions.index:
        if positions.loc[ind, "Strand"] == '+':
            continue
        if positions.loc[ind, "Strand"] == '-':
            tmp = positions.loc[ind, "Three"]
            positions.loc[ind, "Three"] = positions.loc[ind, "Five"]
            positions.loc[ind, "Five"] = tmp
            

    Five_pd_up = pd.DataFrame({"Start": positions.loc[:,"Five"].add(20)}) 
    Five_pd_down = pd.DataFrame({"End": positions.loc[:,"Five"].add(-20)})
    
    Three_pd_up = pd.DataFrame({"Start": positions.loc[:,"Three"].add(20)}) 
    Three_pd_down = pd.DataFrame({"End": positions.loc[:,"Three"].add(-20)})
    
    Five_pd = pd.concat([positions.loc[:,"Chr"], Five_pd_up], axis=1)
    Five_pd = pd.concat([Five_pd, Five_pd_down], axis=1)

    Three_pd = pd.concat([positions.loc[:,"Chr"], Three_pd_up], axis=1)
    Three_pd = pd.concat([Three_pd, Three_pd_down], axis=1)
    
    Five_pd.to_csv(Five_prime_bed, sep="\t")
    Three_pd.to_csv(Three_prime_bed, sep="\t")
    
    

if __name__ == '__main__':
    parser = ArgumentParser(
                        description='Extract 20bp splicing site sequence')
    parser.add_argument('splice_site_file',
                        type=str,
                        help='splice_site_file')
    parser.add_argument('output_bed',
                        type=str,
                        help="output bed file")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

args = parser.parse_args()
extend_splice_sites(args.splice_site_file, args.output_bed)



    