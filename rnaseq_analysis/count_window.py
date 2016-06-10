import pysam
import argparse
import os
import strtime
import multiprocessing
from functools import partial
import sys



def count_window(strand, chromosome, window, step, filenames):
    reads = pysam.AlignmentFile(filenames[0], "rb")
    filename = filenames[1]
    counts = []
    nhits = 0
    names = []
    if chromosome: #for chromosome in reads.header['SQ']:
        index = 0
        while index < chromosome['LN']:
            overlap = []
            name='%s%s:%d-%d' %(chromosome['SN'], strand, index, index + window)
            for read in reads.fetch(chromosome['SN'], index, index + window):
                    #if read.is_read1 and read.is_proper_pair: 
                        #if read.qname not in overlap:
                    
                    overlap.append(read.qname)
                    #elif read.is_proper_pair:
                    #    try:
                    #        mate = reads.mate(read)
                    #        if mate.qname not in overlap:
                    #            overlap.append(mate.qname)
                    #    except:
                    #        pass
            #counts[geneid][(reg[0], reg[1])] = len(set(overlap))
            hits = len(set(overlap))
            nhits += hits
            counts.append([hits, name])
            index += step
    reads.close()
    print "%s: Counted %d hits on strand %s" % (filename, nhits, strand)
    return counts



def count(filesi, fileo, strand, window, step, nproc):
    pool = multiprocessing.Pool(processes=nproc)
    results = [] 
    labels = []
    
    
    reads = pysam.AlignmentFile(filesi[0][0], "rb")
    chromosomes = reads.header['SQ']
    reads.close()
    for chromosome in chromosomes:
        if strand == ".":
            p = partial(count_window, '+', chromosome, window, step)
            results.append(pool.map_async(p, filesi) )
            p = partial(count_window, '-', chromosome, window, step)
            results.append( pool.map_async(p, filesi) )
        else:
            p = partial(count_window, strand, chromosome, window, step)
            results.append(pool.map_async(p, filesi) )
     
    for result_list in results:
        result_list.wait(timeout = 0)
        r = result_list.get()
        for reads in zip(*r):
            name = reads[0][1]
            # check if names are the same
            if sum(x[1] for x in reads if x[1] != name) > 1:
                sys.exit('Bam files mapped on different genomes')
            s = name+"\t"+"\t".join(map(str,[x[0] for x in reads]))
            fileo.write(s+"\n")
    
    pool.close()
    pool.join()
    
def positive_int(val):
    try:
        assert(int(val) >= 0)
    except:
        raise argparse.ArgumentTypeError("'%s' is not a valid positive int" % val)
    return int(val)

def main():
    parser = argparse.ArgumentParser(description='Count reads mapped to genes and exons')
    i = False
    
    parser.add_argument('-b', '--inp', type=str, nargs='+',
                   help='input bam files')
    parser.add_argument('-o', '--out', type=str, 
                   help='output file')
    parser.add_argument('-n', '--nproc', type=positive_int, default=1,
                   help='number of processes')
    parser.add_argument('-l', '--names', type=str, nargs='+',
                   help='input track names', required=False, default=[])

    parser.add_argument("-s", '--strand', help='strand', choices=['+', '-', '.'])
    parser.add_argument("-w", '--window', type=positive_int, required=True)
    parser.add_argument("-t", '--step', type=positive_int, required=True)
    

    args = parser.parse_args()
    
    if args.nproc < 1:
        args.nproc = 1

    if len(args.inp) > len(args.names):
        args.names += args.inp[len(args.names):]

    #if args.e:
    #    print strtime.strtime(), "Counting exon hits on %d slots" % args.nproc
    #elif args.i:
    #    print strtime.strtime(), "Counting gene (with introns) hits on %d slots" % args.nproc
    #elif args.t:
    #    print strtime.strtime(), "Counting transcript hits on %d slots" % args.nproc
    #else:
    #    print strtime.strtime(), "Counting gene (without introns) hits on %d slots" % args.nproc

    fileo = open( args.out, "w+")
    filesi = []
    fileo.write("chrom_strand:start-end")
    for bam in args.inp:
        #try to open before running
        f = open(bam)
        f.close()
        name=args.names.pop(0)
        fileo.write("\t%s"%name)
        if not os.path.exists(bam+".bai"):
            print "Indexing bam file %s"%bam
            pysam.index(bam)
        filesi += [(bam, name)]
        #print bam
    fileo.write("\n")    
    #print "Files: ", filesi 

    count(filesi, fileo, args.strand, args.window, args.step, args.nproc)
    
    #print strtime.strtime()+": Finished"
    
    fileo.close()
    
    
    
if __name__ == "__main__":
   main()
