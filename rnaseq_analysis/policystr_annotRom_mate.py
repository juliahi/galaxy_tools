'''
Created on 28-04-2016
version 0.2 - with annotation of policistronic sides given by Romek Szczesny
@author: irina
'''
import pysam, argparse,sys

def read_annotation(ann, strand):
    annot = {}
    op = open(ann)
    for line in op.readlines():
        lis = line.split('\t')
        if lis[2].split('\n')[0] == strand:
            annot[int(lis[0])]= lis[1]
        else: continue
    return annot

def calc_reads(sam, an, n, genes, b_nr, ind):
#an = {(start,end):'gen_name'}, genes= {gen_name:nr}
    find_mate = {}
    sam_fe = sam.fetch('chrM')
    read_n = 0
    for read in sam_fe:
        read_n +=1 
        name = read.query_name #qname
        
        #if read_n % 10000 ==0: sys.stderr.write("%d\n"%read_n)
        
        if read.mate_is_unmapped:
            pos = read.get_reference_positions()
            position = []
            for key in an.keys():
                if key - n in pos and key + n in pos:
                    #print "Dodaje", an[key], name
                    position.append(an[key])
                else: pass
            if len(position) == 0: 
                continue
            po_name = '_'.join(position)
            if po_name in genes:
                genes[po_name][ind] += 1
            else: 
                genes[po_name] = [0]*b_nr
                genes[po_name][ind] += 1 
            
        elif name in find_mate:
            pos_mate = read.get_reference_positions() # para znaleziona
            pos = find_mate[name].get_reference_positions()
            position = []
            #if len(position)> 0: print "kasuje", len(position)
            for key in an:
                #print pos[0]+n,int(key), pos[-1]-n, "|", pos_mate[0]+n, int(key), pos_mate[-1]-n
                #if pos[0]+n <= int(key) <=  pos[-1]-n or pos_mate[0]+n <= int(key) <=  pos_mate[-1]-n:
                if (key - n in pos and key + n in pos) or (key - n in pos_mate and key + n in pos_mate):
                    #print "Dodaje", an[key], name
                    position.append(an[key])
                else: pass
            #if len(position)>4: print position, "MATE", name, read.get_reference_positions(), "sec", find_mate[name].get_reference_positions()
            if len(position) == 0: 
                del find_mate[name]
                continue
            po_name = '_'.join(position)
            if po_name in genes:
                genes[po_name][ind] += 1
            else: 
                genes[po_name] = [0]*b_nr
                genes[po_name][ind] += 1 
            del find_mate[name]
        else: find_mate[name] = read
    #print "Zostalo", len(find_mate.keys()), read_n
    
    for name in find_mate:
            pos = find_mate[name].get_reference_positions()
            position = []
            for key in an:
                if key - n in pos and key + n in pos:
                    #print "Dodaje", an[key], name
                    position.append(an[key])
                else: pass
            #if len(position)>4: print position, "MATE", name, read.get_reference_positions(), "sec", find_mate[name].get_reference_positions()
            if len(position) == 0:
                continue
            po_name = '_'.join(position)
            if po_name in genes:
                genes[po_name][ind] += 1
            else: 
                genes[po_name] = [0]*b_nr
                genes[po_name][ind] += 1 
            
    #print len(find_mate.keys()), read_n   
                  
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
    #print annot_d
    N = args.Overlap
    ge = {}
    bam_n = len(args.BamFile)
    i = 0
    
    #print 'Annot', annot_d, 'Overl', N, 'genes', ge, 'File_nr', bam_n, 'now_nr', i
    for bam_file in args.BamFile: 
        samfile = pysam.AlignmentFile(bam_file, "rb")
        
        sam_fe = samfile.fetch('chrM')
        read_n = 0
        for read in sam_fe:
            read_n +=1
        print "READS = ", read_n
        
        ge = calc_reads(samfile, annot_d, N, ge, bam_n, i)
        i += 1
        
    if len(args.Column) == 0:
       header = [' '] + args.BamFile
    else: header = [' '] + args.Column
       
    if args.Out == "":
       print '\t'.join(header)
       #print sorted(ge.keys())
       for gen_k in sorted(ge.keys()):
           li = '\t'.join(map(str,ge[gen_k]))
           print gen_k+'\t'+li
    else: 
       outfile = open(args.Out, 'w')
       line = '\t'.join(header)
       outfile.write(line + '\n')
       for gen_k in sorted(ge.keys()):
           li = '\t'.join(map(str,ge[gen_k]))
           outfile.write(gen_k+'\t'+li+'\n')
 
    
      
        
      
