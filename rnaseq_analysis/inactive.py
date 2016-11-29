#!/usr/bin/python

import pysam
import argparse

import numpy as np
import pandas

import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns

import sys
from os import path
sys.path.insert(0,  path.dirname( path.dirname( path.abspath(__file__) ) ) )

import count_for_deseq_par_list as countlib

from math import log

def plot_hist(args, counts, d):
    
    #plot
    n = len(counts.columns)
    chrt = 0
    fig = plt.figure()
    
    q = counts.quantile(q=0.99, axis=0)
    qs = counts.quantile(q=0.9, axis=0)
    #q1 = max(q1, )
    bound = float(max(q))
    #bound = int(counts.values.max()+1)
    #if bound > 0:
    #counts = counts.applymap(lambda x: min(x, bound))
    
    for column in counts:
        chrt += 1 
        ax = fig.add_subplot( 2,(n+1)/2, chrt)
        ax.set(title=column)    
        col = list(counts[column])
        col = [log(x+1) for x in col]
        q1 = float(qs[column])
        #print list( np.linspace(0,q1,10, endpoint=False) ) 
        #print list( np.linspace(q1, bound+1, 50) )
        bins = list( np.linspace(0,q1,50, endpoint=False) ) + list( np.linspace(q1, bound+1, 10) )
        bins = list( np.linspace(0, bound, 100))
        dpl = sns.distplot(col, label=column, rug=False, norm_hist=False, bins=bins, ax=ax, kde_kws={'bw': 3}, kde=False)
        dpl.set_xlim([0, bound])
        dpl.set(xlabel="log(TPM)", ylabel="Number of genes")
        #break
    
    plt.savefig(args.output[2], format='pdf')
    

def main():
    parser = argparse.ArgumentParser(description='Find inactive genes')
    
    parser.add_argument('-i', '--input', type=str, required=True,
                   help='input deseq result')
    parser.add_argument('-a', '--annotation', type=str, required=True,
                   help='input annotation')
    parser.add_argument('-t', '--counttype', type=str, choices=['e', 't', 'g', 'i'], required=True,
                   help='what to count: e=exons t=transcripts g=genes i=genes with introns')
    parser.add_argument('-o', '--output', type=str, nargs=3, required=True,
                   help='output files: inactive genes, active genes counts, pdf with histograms')
    
    args = parser.parse_args()
    
    
    

    genesp, namesp = countlib.gen_genes(args.annotation, args.counttype, '+')

    genesm, namesm = countlib.gen_genes(args.annotation, args.counttype, '-')
    d = dict(zip(namesp+namesm, [g.length(args.counttype) for g in genesp+genesm]))

    counts  = pandas.read_csv(args.input, sep='\t', index_col=0)
    #change counts to TPMs! first normalize by length, then by mapped transcripts
    for column in counts:
        col = list(counts[column])
        
        lengths = [d[x] for x in list(counts.index) if x in d]
        counts[column] = [float(x)/y for x, y in zip(col, lengths)]
    for column in counts:
        col = list(counts[column])
        colsum = sum(col)
        counts[column].apply(lambda x: x*(10**6)/colsum + 1 )


    inputf = open(args.input)
    header = inputf.readline()
    
    with open(args.output[0], 'w+') as inactive:
        with open(args.output[1], 'w+') as active:
            active.write(header)
            for line in inputf:
               fields = line.strip().split('\t')
               if fields[0] in d:
                   if sum([1 for i in range(len(fields)-1) if list(counts.loc[fields[0]])[i] >= 2]) > 0: #FPKM > 2 -> active
                        active.write(line)
                   else:
                        inactive.write(fields[0]+'\n')
               else:
                   print fields[0] + ' not found in annotation'
    inputf.close()
    
    plot_hist(args, counts, d)


             
    
if __name__ == "__main__":
   main()


