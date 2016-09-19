#!/usr/bin/python

import pysam
import argparse
import strtime
import os

def count_all(bamfile):
    #all reads: mapped + unmapped
    #return reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfile) ])
    
    #count mapped reads
    return reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(bamfile) ])

def coverage(inputs, output, normalize):
    inputfiles = []
    for inp in inputs:
        #if not os.path.exists(inp+".bai"):
        #    print "Indexing %s" % inp
        #    pysam.index(inp)
        inputfiles.append( pysam.AlignmentFile(inp, "rb") )
    outputfile = open(output, "w")
    
    steps = []
    for inp, number in zip(inputs, normalize):
            print number
            steps.append( 1000000.0 / number )
    print  steps
    
    
    
    for seqs in inputfiles[0].header['SQ']:
        chromosome = seqs['SN']
        length = seqs['LN']
        #print "Counting chromosome", chromosome
        
        starts = []
        ends = []
        for inp, step in zip(inputfiles, steps):
            for read in inp.fetch(chromosome):
                for block in read.get_blocks():
                    starts.append((block[0], step))
                    ends.append((block[1], step))
                    
        # mozna zmienic na efektywne dodawanie do listy
        starts.sort()
        ends.sort()
        x = 0
        y = 0
        size = 0
        pocz = 0
        kon = 0
        oldsize = 0
        lines = []
        while True:
            while x < len(starts) and starts[x][0] == kon:
                size += starts[x][1]
                x += 1
                #print x,y, starts[x], ends[y]
            while y < len(ends) and ends[y][0] == kon:
                size -= ends[y][1]
                y += 1
                #print x,y, starts[x], ends[y]
            if x < len(starts) and starts[x][0] < ends[y][0]:
                kon = starts[x][0]
                #size += starts[x][1]
                #x += 1
                
            elif y < len(ends):
                kon = ends[y][0]
                #size -= ends[y][1]
                #y += 1
            else:
                kon = length
           
            lines.append( [ pocz, kon, size] )
            pocz = kon
            if kon == length:
                break
            
        i = 0
        while i < len(lines):
            if i < len(lines)-1 and lines[i][2]==lines[i+1][2]:
                lines[i+1][0] = lines[i][0]
            else:
                outputfile.write("%s\t%d\t%d\t%f\n" % (chromosome, lines[i][0], lines[i][1], lines[i][2]) )
            i += 1
    outputfile.close()

def main():
    parser = argparse.ArgumentParser(description='Coverage BedGraph')
    
    parser.add_argument('-i', '--input', type=str, nargs='+', required=True,
                   help='input bam files (plus_reads minus_reads respectively)')
    parser.add_argument('-o', '--output', type=str, nargs=2, required=True,
                   help='output file')
    
    #parser.add_argument("-s", '--strand', type=str, help='strand')
    parser.add_argument('-n', type=int, nargs='+',
                   help='normalize by this number')
    
    args = parser.parse_args()
    
    if len(args.input) != 2*len(args.n):
        parser.error("Wrong number of values to normalize by:%d files, %d values. Should be one for each pair")
    
    print strtime.strtime() 
      
    print "Counting coverage of files:", args.input

    open(args.output[0], "w").close()
    open(args.output[1], "w").close()
    for i in args.input:
        open(i).close()
    coverage(args.input[::2], args.output[0], args.n)
    
    coverage(args.input[1::2], args.output[1], args.n)
    
    print strtime.strtime()+": Finished"

    
    
if __name__ == "__main__":
   main()


