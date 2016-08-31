#!/usr/bin/python

import pysam
import argparse


def main():
    parser = argparse.ArgumentParser(description='Add gene names to DESeq output')
    
    parser.add_argument('-i', '--input', type=str, required=True,
                   help='input deseq result')
    parser.add_argument('-o', '--output', type=str, required=True,
                   help='output file')
    
    args = parser.parse_args()
    
    with open(args.input) as f:
        with open(args.output, 'w') as output_file:
            header = f.readline()
            output_file.write(header.strip() + "\tgeneName\n" )
            for line in f:
                fields = line.strip().split('\t')
                if len(fields[0].split('_')) > 2:
                    name = fields[0].split('_')[1:-1]
                    name = "_".join(name)
                else:
                    name = fields[0]
                fields += [name]
                output_file.write( "\t".join(fields)+'\n')
    
    
if __name__ == "__main__":
   main()


