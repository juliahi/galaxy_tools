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
    if normalize:
        #count coverage in reads/1M reads
        for inp in inputs:
            print count_all(inp)
            steps.append( 1000000.0 / count_all(inp) )
    else:
        steps = [1] * len(inputs)
    print "Mean coverages:", steps
    
    
    
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
                   help='input bam files')
    parser.add_argument('-o', '--output', type=str, required=True,
                   help='output file')
    
    #parser.add_argument("-s", '--strand', type=str, help='strand')
    parser.add_argument('-n', action="store_true", 
                   help='normalize by read count in bam file')
    
    args = parser.parse_args()
    
    print strtime.strtime() 
    if args.n:
        print "Counting normalized coverage of files:", args.input
    else:        
        print "Counting coverage of files:", args.input

    open(args.output, "w").close()
    for i in args.input:
        open(i).close()
    coverage(args.input, args.output, args.n)
    
    print strtime.strtime()+": Finished"

    
    
if __name__ == "__main__":
   main()


