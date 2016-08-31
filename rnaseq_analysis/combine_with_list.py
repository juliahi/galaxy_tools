#!/usr/bin/python

import pysam
import argparse


def main():
    parser = argparse.ArgumentParser(description='Add gene names to DESeq output')
    
    parser.add_argument('-i', '--input', type=str, required=True,
                   help='input file to extend')
    parser.add_argument('-l', '--lists', type=str, nargs='+', 
                   required=True, help='files with values')
    parser.add_argument('-c', '--fields', type=int, nargs='+', required=True,
                   help='numbers of fields to compare')
    parser.add_argument('-n', '--names', type=str, nargs='+', required=True,
                   help='names of fields to compare')
    parser.add_argument('-o', '--output', type=str, required=True,
                   help='output file')
    
    
    args = parser.parse_args()
    
    if len(args.fields) != len(args.lists):
        parser.error("Number of given lists of values and fields should be equal")
    
    
    with open(args.input) as f:
        value_lists = []
        for filename in args.lists:
            with open(filename) as f2:
                value_lists.append([l.strip() for l in f2])
            
        with open(args.output, 'w') as output_file:
            header = f.readline()
            output_file.write(header.strip() + "\t" + "\t".join(args.names)+"\n" )
            for line in f:
                fields = line.strip().split('\t')
                for values, colnum in zip(value_lists, args.fields):
                    if len(fields) >= colnum:
                        if fields[colnum-1] in values:
                            fields.append('+')
                        else:
                            fields.append('-')
                output_file.write( "\t".join(fields)+'\n')
    
    
if __name__ == "__main__":
   main()


