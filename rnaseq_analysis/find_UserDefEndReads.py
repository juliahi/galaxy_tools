'''
Created on 19-11-2016
version 0.1 - extract reads with 3' end defined by the user
@author: irina
'''
import pysam, argparse, sys

def find_end(sam, out, end):
    samf_end = pysam.AlignmentFile(out, "wb", template=sam)
    le = len(end)
    print end
    sam_fe = sam.fetch()
    for read in sam_fe:
        if end ==  read.query_sequence[-le:]:
            #print read.query_sequence
            samf_end.write(read)
        else: pass
            

if __name__ == '__main__':
    
    optparser = argparse.ArgumentParser()
    optparser.add_argument('-f', type = str, dest="BamFile", default = "", help = "Bam file")
    optparser.add_argument('-e', type = str, dest="SeqEnd", default = "", help = "Sequence that we want to be at the 3' end of read. Only capital letters without space between them: CCACCA")
    optparser.add_argument('-o', type = str, default = "reads_with_userDefEnd.bam", dest="Output", help = "Output filename")
    
    args = optparser.parse_args()
    if len(sys.argv) ==1:
        print optparser.print_help()
        sys.exit(1)
    
    if args.SeqEnd == "":
        print "Please, give the sequence you are looking for. Use -e option!!!!"
        #print optparser.print_help()
        sys.exit(1)
        
    samfile = pysam.AlignmentFile(args.BamFile, "rb")
    find_end(samfile, args.Output, args.SeqEnd)
