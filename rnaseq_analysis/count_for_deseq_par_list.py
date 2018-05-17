import pysam
import argparse
import os
#import strtime
import multiprocessing
from functools import partial
import sys

class Gene:
    name = ""
    strand = "."
    chrom = ""
    start = 0
    end = 0
    exons = set()
    def __init__(self, name, chrom, strand, start, end, exons):
        self.name = name
        self.strand = strand
        self.chrom = chrom
        self.start = start
        self.end = end
        self.exons = set(exons)
    def update(self, exons):
        self.exons.update(exons)
    def length(self, typ):
        if typ == 'g':
            return self.end-self.start
        return sum([e-s for e,s in exons])
    

def gen_genes(annotation, counttype, strand):
    """ strand shoult be + or -. Doesnt work for both strands"""
    """counttype: t for transcripts, i for genes with introns included, 
    g for genes without introns, e for exons"""
    genes = []
    
    names = []
    gene_ids = []
    prev_chrom = None
    start_pos = 0

    id_count = 0
    for line in open(annotation, "r"):
        if line[0] == "#": continue
        g = line.strip().split('\t')
        if g[2] != strand: continue

        exons = zip(map(int, [x for x in g[5].split(",") if x != '']), map(int, [x for x in g[6].split(",") if x != '']))
        start = int(g[3])
        end = int(g[4])
        if counttype=="t":
            genes += [Gene(g[0]+'_'+g[7], g[1], g[2], start, end, exons)]
            #names.append(g[0]+"_"+g[7])
            continue
        if prev_chrom != g[1]:
            start_pos += len(gene_ids)
            gene_ids = []
            prev_chrom = g[1]
        
        if len(g) <= 8:
            g.append(g[7]+'_'+g[1])
        pos = None
        if gene_ids != [] and g[8] == gene_ids[-1]:
            pos = len(gene_ids)-1
        elif g[8] in gene_ids:
            pos = gene_ids.index(g[8])
        if pos is not None:
            gene = genes[start_pos + pos]
            print gene.name, pos, g[8]
            gene.start=min(gene.start, start)
            gene.end=max(gene.end, end)
            if counttype != "i":
                     gene.update(exons)
        else:
            genes += [Gene(g[8]+'_'+g[7], g[1], g[2], start, end, exons)]
            gene_ids.append( g[8] )
            id_count += 1
        
            
    if counttype == "g" or counttype == "i" or counttype == "t":
        for gene in genes:
            if counttype == "g":
                ex = sorted(list(gene.exons))
                regions = [list(ex.pop(0))]
                for reg in ex:
                    if reg[0] <= regions[-1][1]:
                        if reg[1] > regions[-1][1]:
                            regions[-1][1] = reg[1]
                    else:
                        regions.append(list(reg))
                gene.exons = regions
            name="%s_%s%s:%d-%d"%(gene.name, gene.chrom, gene.strand, gene.start, gene.end)
            names.append(name) 
    if counttype == "e":
        exons = []
        for gene in genes:
            for ex in sorted(list(gene.exons)):
                name="%s_%s%s:%d-%d"%(gene.name, gene.chrom, gene.strand, ex[0], ex[1])
                names.append(name) 
                exons.append(Gene(gene.name, gene.chrom, gene.strand, ex[0], ex[1], [ex]))
        return (exons, names)
    return (genes, names)



def count_exons(genes, filenames):
    reads = pysam.AlignmentFile(filenames[0], "rb")
    filename = filenames[1]
    counts = []
    ngenes = 0
    nhits = 0
    
    for gene in genes:
        #counts[geneid] = {}
        regions = sorted(list(gene.exons))
                
        
        for reg in regions:
                #geneid="%s_%s%s:%d-%d"%(gene.name, gene.chrom, gene.strand, reg[0], reg[1])
                if gene.chrom not in reads.references:
                    print "Warning: There is no chromosome %s in BAM file reference! " % gene.chrom
                    hits = 0
                else:    
                    overlap = []
                    for read in reads.fetch(gene.chrom, reg[0], reg[1]):
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
                counts.append(hits)
                ngenes += 1      
                
    reads.close()
    if ngenes == 0:
        sys.exit("No genes in file %s"%filename)
    else:
        print "%s: Counted %d hits for %d exons on strand %s" % (filename, nhits, ngenes, genes[0].strand)
    return counts

def count_genes(genes, counttype, filenames):
    reads = pysam.AlignmentFile(filenames[0], "rb")
    filename = filenames[1]
    counts = []
    ngenes = 0
    nhits = 0
    
    for gene in genes:
        #geneid="%s_%s%s:%d-%d"%(gene.name, gene.chrom, gene.strand, gene.start, gene.end)
        if counttype=="i":
            regions = [[gene.start, gene.end]]
        else:
            regions = gene.exons
        overlap = []
        for reg in regions:
                if gene.chrom not in reads.references:
                    print "Warning: There is no chromosome %s in BAM file reference! " % gene.chrom
                else:    
                    for read in reads.fetch(gene.chrom, reg[0], reg[1]):
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
                #counts[geneid] = len(set(overlap))
        hits=len(set(overlap))
        nhits += hits
        counts.append(hits)
        ngenes += 1
    reads.close()
    if ngenes == 0:
        print "No genes in file %s"%filename
    elif counttype == "g":
        print "%s: Counted %d hits for %d genes on strand %s" % (filename, nhits, ngenes, genes[0].strand)
    elif counttype == "i":
        print "%s: Counted %d hits for %d genes (with introns) on strand %s" % (filename, nhits, ngenes, genes[0].strand)
    elif counttype == "t":
        print "%s: Counted %d hits for %d transcripts on strand %s" % (filename, nhits, ngenes, genes[0].strand)
    elif counttype == "e":
        print "%s: Counted %d hits for %d exons on strand %s" % (filename, nhits, ngenes, genes[0].strand)
    
    return counts
    #for gid, n in counts.iteritems():
    #    if not count_exons:
    #        outputfile.write("%s\t%d\n"%(gid, n))
    #        ngenes += 1
    #        nhits += n
    #    else:
    #        for ex, c in n.iteritems():
    #            outputfile.write("%s_%d-%d\t%d\n"%(gid, ex[0], ex[1], c))
    #            nhits += c
    #            ngenes += 1
   

def count(annotation, filesi, fileo, strand, counttype, nproc):
    if strand == ".":
        strands = ["+", "-"]
    else:
        strands = [strand]
    
    p = partial(gen_genes, annotation, counttype)
    #genesl = pool.apply_async(p, strands)
    #genesl.wait()
    genesl = map(lambda(x): gen_genes(annotation, counttype, x), strands)
    

    if len(genesl) == 0 or len(genesl[0][0]) == 0:
        return 
    pool = multiprocessing.Pool(processes=nproc)
    results = [] 
    labels = []
    for genes, names in genesl:
        if counttype == "e":
            results.append( pool.map_async(partial(count_exons, genes), filesi) )
        else:
            #results = pool.apply_async(lambda(x): count_genes(genes, x[0], counttype, x[1]), filesi)
            p = partial(count_genes, genes, counttype)
            results.append( pool.map_async(p, filesi) )
        labels.append(names)

    for result, names in zip(results, labels):
        result.wait()
        r = result.get()
        for reads in zip(*[names]+r):
            s = "\t".join(map(str,reads))
            fileo.write(s+"\n")
    
    pool.close()
    pool.join()

def main():
    parser = argparse.ArgumentParser(description='Count reads mapped to genes and exons')
    i = False
    
    parser.add_argument('annotation', type=str, 
                   help='annotation file')
    parser.add_argument('-b', '--inp', type=str, nargs='+',
                   help='input bam files')
    parser.add_argument('-o', '--out', type=str, 
                   help='output file')
    parser.add_argument('-n', '--nproc', type=int, default=1, 
                   help='number of processes')
    parser.add_argument('-l', '--names', type=str, nargs='+',
                   help='input track names', required=False, default=[])

    parser.add_argument("-s", '--strand', help='strand', choices=['+', '-'])
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-i', action="store_true", 
                   help='including introns in gene hits count')
    group.add_argument('-e', action="store_true", 
                   help='count exon hits')
    group.add_argument('-t', action="store_true", 
                   help='count transcript hits')
    
    args = parser.parse_args()
     
    if args.nproc < 1:
        args.nproc = 1

    if len(args.inp) > len(args.names):
        args.names += args.inp[len(args.names):]
    print args.names
    #if args.e:
    #    print strtime.strtime(), "Counting exon hits on %d slots" % args.nproc
    #elif args.i:
    #    print strtime.strtime(), "Counting gene (with introns) hits on %d slots" % args.nproc
    #elif args.t:
    #    print strtime.strtime(), "Counting transcript hits on %d slots" % args.nproc
    #else:
    #    print strtime.strtime(), "Counting gene (without introns) hits on %d slots" % args.nproc

    filea = open(args.annotation, "r")
    filea.close()
    fileo = open( args.out, "w+")
    filesi = []
    #fileo.write("#")
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

    #if args.e:
    #    count_exons(filea, filesi, fileo, args.strand)
    #else:
    if True:
        if args.i:
            counttype = "i" #genes with only intron hits
        elif args.t:
            counttype = "t" #transcripts
        elif args.e:
            counttype = "e" #exons
        else:
            counttype = "g" #genes
        count(args.annotation, filesi, fileo, args.strand, counttype, args.nproc)
        #count_genes(filea, filesi, fileo, args.strand, counttype)
    
    #print strtime.strtime()+": Finished"
    
    #filea.close()
    #for f in filesi:
    #    f.close()
    #fileo.close()
    
    
    
if __name__ == "__main__":
   main()
