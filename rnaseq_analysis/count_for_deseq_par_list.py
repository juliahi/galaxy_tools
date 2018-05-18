import pysam
import argparse
import os
#import strtime
import multiprocessing
from functools import partial
import sys

class Gene:
    # name = ""
    # strand = "."
    # chrom = ""
    # start = 0
    # end = 0
    # exons = set()
    def __init__(self, name, chrom, strand, start, end, exons):
        self.name = name
        self.strand = strand
        self.chrom = chrom
        self.start = start      # 0-based!
        self.end = end          # 0-based included?!
        self.exons = set(exons)
    def update(self, exons):
        self.exons.update(exons)
    def length(self, typ):
        if typ == 'g':
            return self.end-self.start
        return sum([e-s for e,s in exons])

    def join_exons(self):
        ex = sorted(list(self.exons))
        regions = [list(ex.pop(0))]
        for reg in ex:
            if reg[0] <= regions[-1][1]:
                if reg[1] > regions[-1][1]:
                    regions[-1][1] = reg[1]
            else:
                regions.append(list(reg))
        self.exons = [tuple(x) for x in regions]

    def get_introns(self):
        introns = []
        for ex1, ex2 in zip(self.exons[:-1], self.exons[1:]):
            introns.append((ex1[1], ex2[0]))
        return introns


def read_line(line, strand):
    if line[0] == "#": return None
    g = line.strip().split('\t')
    if g[2] != strand: return None

    tr_id = g[0]
    chromosome = g[1]
    exons = zip(map(int, [x for x in g[5].split(",") if x != '']), map(int, [x for x in g[6].split(",") if x != '']))
    start = int(g[3])
    end = int(g[4])

    gene_name=g[7]
    if len(g) <= 8:
        gene_id = g[7] + '_' + chromosome
    else: gene_id = g[8]

    return tr_id, chromosome, strand, start, end, exons, gene_name, gene_id


# Generate list of Gene objects for annotation, for specific counting type
def read_annotation(annotation, counttype, strand):
    """ strand should be + or -. Doesnt work for both strands"""
    """ counttype[0] can be 'genes', 'transcripts' or 'exons' or 'introns 
        -- returns one Gene entry for each gene/transcript/gene exon/gene intron
    counttype[1] can be 'all' or 'exons' or 'introns' 
        -- for genes/transcripts means what to count as 
    """
    genes = []

    gene_ids = []
    prev_chrom = None
    start_pos = 0

    for line in open(annotation, "r"):
        line_result = read_line(line, strand)
        if not line_result: continue
        else:
            tr_id, chromosome, strand, start, end, exons, gene_name, gene_id = line_result

        if counttype[0]=="transcripts":
            if counttype[1] == 'all': exons = [(start, end)]  # TODO: check this!
            genes += [Gene(tr_id+'_'+gene_name, chromosome, strand, start, end, exons)]
            #names.append(g[0]+"_"+g[7])
            continue
        else:
            if prev_chrom != chromosome:
                start_pos += len(gene_ids)
                gene_ids = []
                prev_chrom = chromosome

            # check if variant of this gene was already loaded
            pos = None
            if gene_ids != [] and gene_id == gene_ids[-1]:
                pos = len(gene_ids)-1
            elif gene_id in gene_ids:
                pos = gene_ids.index(gene_id)

            if pos is not None:     # update gene
                gene = genes[start_pos + pos]
                # print gene.name, pos, gene_id
                gene.start=min(gene.start, start)
                gene.end=max(gene.end, end)
                if counttype != ("genes", 'all'):
                    gene.update(exons)
                else:
                    gene.exons = [(gene.start, gene.end)]
            else:   # add gene
                if counttype == ("genes", 'all'): exons = [(start, end)]   # TODO: check this!
                genes += [Gene(gene_id+'_'+gene_name, chromosome, strand, start, end, exons)]
                gene_ids.append(gene_id )

    names = []

    if counttype[0] == "genes" or counttype[0] == "transcripts":
        for gene in genes:
            if counttype[0] == "genes": gene.join_exons()
            if counttype[1] == 'all': gene.exons = [(gene.start, gene.end)]
            if counttype[1] == 'introns': gene.exons = gene.get_introns()
            name = "%s_%s%s:%d-%d" % (gene.name, gene.chrom, gene.strand, gene.start, gene.end)
            names.append(name)
        return genes, names

    if counttype[0] == "exons":
        exons = []
        for gene in genes:
            gene.join_exons()
            for ex in sorted(list(gene.exons)):
                name="%s_%s%s:%d-%d"%(gene.name, gene.chrom, gene.strand, ex[0], ex[1])
                names.append(name)
                exons.append(Gene(gene.name, gene.chrom, gene.strand, ex[0], ex[1], [ex]))
        return exons, names

    if counttype[0] == "introns":
        introns = []
        for gene in genes:
            gene.join_exons()
            for ex in gene.get_introns():
                name="%s_%s%s:%d-%d"%(gene.name, gene.chrom, gene.strand, ex[0], ex[1])
                names.append(name)
                introns.append(Gene(gene.name, gene.chrom, gene.strand, ex[0], ex[1], [ex]))
        return introns, names


def count_genes(genes, counttype, filenames):
    reads = pysam.AlignmentFile(filenames[0], "rb")
    filename = filenames[1]
    counts = []
    ngenes = 0
    nhits = 0
    
    for gene in genes:
        overlap = []
        for reg in gene.exons:
            if gene.chrom not in reads.references:
                print "Warning: There is no chromosome %s in BAM file reference! " % gene.chrom
            else:
                for read in reads.fetch(gene.chrom, reg[0], reg[1]):
                    #if read.is_read1 and read.is_proper_pair:
                        #if read.qname not in overlap:

                    # reads are 0-based, but use reg[0]-1 to get overlaps with ends of reads
                    if read.get_overlap(reg[0], reg[1]) > 0:        # CRITICAL BUG FIX!!!
                        overlap.append(read.qname)
        hits=len(set(overlap))
        nhits += hits
        counts.append(hits)
        ngenes += 1
    reads.close()
    if ngenes == 0:
        print "No genes in file %s"%filename
    else:
        print "%s: Counted %d hits for %d %s (%s) on strand %s" % (filename, nhits, ngenes,
                                                                   counttype[0], counttype[1], genes[0].strand)
    return counts


# async counting
def count(annotation, filesi, fileo, strand, counttype, nproc):
    if strand == ".": strands = ["+", "-"]
    else: strands = [strand]


    annotation_list = map(lambda(x): read_annotation(annotation, counttype, x), strands)

    if len(annotation_list) == 0 or len(annotation_list[0][0]) == 0:
        return 
    pool = multiprocessing.Pool(processes=nproc)
    results = [] 
    labels = []
    for genes, names in annotation_list:
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

    parser.add_argument('annotation', type=str, help='annotation file')
    parser.add_argument('-b', '--inp', type=str, nargs='+', help='input bam files')
    parser.add_argument('-o', '--out', type=str, help='output file')
    parser.add_argument('-n', '--nproc', type=int, default=1, help='number of processes')
    parser.add_argument('-l', '--names', type=str, nargs='+',
                   help='input track names', required=False, default=[])
    parser.add_argument("-s", '--strand', help='strand', choices=['+', '-'])

    # group1 = parser.add_mutually_exclusive_group()
    # group1.add_argument('-e', action="store_true",
    #                help='count exon/intron hits')
    # group1.add_argument('-t', action="store_true",
    #                help='count transcript hits')
    # group1.add_argument('-g', action="store_true",
    #                help='count gene hits')
    #
    # group2 = parser.add_mutually_exclusive_group()
    # group2.add_argument('-i', action="store_true",
    #                help='only introns')
    # group2.add_argument('-e', action="store_true",
    #                help='only exons')
    # group2.add_argument('-t', action="store_true",
    #                help='all (begin-end)')

    parser.add_argument('-f', '--feature', choices=('genes', 'transcripts', 'exons', 'introns'), required=True,
                        help='Output one value per gene/transcript/exon/intron')
    parser.add_argument('-i', '--include', choices=('all', 'exons', 'introns'), default='exons',
                        help='For genes count only reads mapped to exons/introns/both.  \
                              Has no efect if feature==exon or intron ')

    args = parser.parse_args()
     
    if args.nproc < 1: args.nproc = 1

    if len(args.inp) > len(args.names):
        args.names += args.inp[len(args.names):]
    print "Counting samples named:", args.names

    # try to open before running
    filea = open(args.annotation, "r")
    filea.close()
    fileo = open( args.out, "w+")
    filesi = []
    for bam in args.inp:
        f = open(bam)
        f.close()
        name=args.names.pop(0)
        fileo.write("\t%s"%name)
        if not os.path.exists(bam+".bai"):
            print "Indexing bam file %s"%bam
            pysam.index(bam)
        filesi += [(bam, name)]
    fileo.write("\n")

    counttype=(args.feature, args.include)
    count(args.annotation, filesi, fileo, args.strand, counttype, args.nproc)

    fileo.close()

if __name__ == "__main__":
   main()
