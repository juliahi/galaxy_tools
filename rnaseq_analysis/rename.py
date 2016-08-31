#!/usr/bin/python

import pysam
import argparse


def main():
    parser = argparse.ArgumentParser(description='Add gene names to DESeq output')
    
    parser.add_argument('-i', '--input', type=str, required=True,
                   help='input file to modify')
    parser.add_argument('-d', '--dict', type=str, required=True,
                   help='dictionary file')
    parser.add_argument('-c', '--columns', type=int, nargs=2, 
                   required=True, help='column numbers')
    parser.add_argument('-o', '--output', type=str, required=True,
                   help='output file')
    
    
    args = parser.parse_args()
    
    
    with open(args.input) as f:
        values = {}
        with open(args.dict) as f2:
           for line in f2:
               fields = line.strip().split('\t')
               if len(fields) >= args.columns[0] and len(fields) >= args.columns[1]:
                   values[fields[args.columns[0]-1]] = fields[args.columns[1]-1]
        with open(args.output, 'w') as output_file:
            for line in f:
                name = line.strip()
                name2 = name.split('.')[0].split('_')[0]
                if name in values:
                    output_file.write( values[name] + '\n')
                elif name2 in values:
                    output_file.write( values[name2] + '\n')
                else:
                    output_file.write( line )
    
if __name__ == "__main__":
   main()


