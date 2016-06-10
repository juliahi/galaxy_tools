'''
Created on 28-04-2016

@author: irina
'''
import pysam, argparse,sys

def read_annotation(ann, strand):
    annot = {}
    op = open(ann)
    for line in op.readlines():
        lis = line.split('\t')
        if lis[2] == strand and lis[1] == 'chrM':
            annot[(lis[3], lis[4])]= lis[8].split('\n')[0]+'_'+lis[7] + '_' + lis[1] + lis[2] + ':' + lis[3] + '-' +lis[4]
        else: continue
    return annot

def calc_reads(sam, an, n, genes, b_nr, ind):
#an = {(start,end):'gen_name'}, genes= {gen_name:nr}
    find_mate = {}
    sam_fe = sam.fetch('chrM')
    for read in sam_fe:
        pos = read.get_reference_positions()
        name = read.query_name #qname
        for key in an.keys():
            if (int(key[0]) <= pos[0] < int(key[-1])-n and pos[-1] > int(key[-1])) or (int(key[0])+n < pos[-1] <= int(key[-1]) and pos[0] < int(key[0])): 
              #print key, pos
              if an[key] in find_mate:
                  if name in find_mate[an[key]]: continue # exclude mate that cross start or end of gene annotation
                  else: find_mate[an[key]].append(name)
              else: find_mate[an[key]] = [name]
              if an[key] in genes:
                  genes[an[key]][ind] += 1
              else: 
                  genes[an[key]] = [0]*b_nr
                  genes[an[key]][ind] += 1
                  
    return genes

if __name__ == '__main__':
    
    optparser = argparse.ArgumentParser()
    optparser.add_argument('-f', type = str, dest="BamFile", help = "Bam file/s", nargs = '+')
    optparser.add_argument('-a', type = str, dest="Annot", help = "Annotation file in txt format")
    optparser.add_argument('-s', type = str, dest="Strand", help = "Strand + or -", choices = ['+', '-'])
    optparser.add_argument('-o', type = str, default = "", dest="Out", help = "Output name",  required=False)
    optparser.add_argument('-c', type = str, default = [], dest="Column", help = "Input track names", nargs = '+', required=False)
    optparser.add_argument('-n', type = int, default = 5, dest="Overlap", help = "Reads nucleotides that should overlap the gene border",  required=False)
    
    args = optparser.parse_args()
    if len(sys.argv) ==1:
        print optparser.print_help()
        sys.exit(1)
    
    #print opts.Annot
    annot_d = read_annotation(args.Annot, args.Strand)
    N = args.Overlap
    ge = {}
    bam_n = len(args.BamFile)
    i = 0
    
    for bam_file in args.BamFile: 
        samfile = pysam.AlignmentFile(bam_file, "rb")
        ge = calc_reads(samfile, annot_d, N, ge, bam_n, i)
        i += 1
        
    if len(args.Column) == 0:
       header = [' '] + args.BamFile
    else: header = [' '] + args.Column
       
    if args.Out == "":
       print '\t'.join(header)
       for gen_k in sorted(ge.keys()):
           li = '\t'.join(map(str,ge[gen_k]))
           print gen_k+'\t'+li+'\n'
    else: 
       outfile = open(args.Out, 'w')
       line = '\t'.join(header)
       outfile.write(line + '\n')
       for gen_k in sorted(ge.keys()):
           li = '\t'.join(map(str,ge[gen_k]))
           outfile.write(gen_k+'\t'+li+'\n')
 
    
      
        
      
