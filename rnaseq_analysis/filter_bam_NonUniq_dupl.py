'''
Created on 13-09-2016
version 0.1 - extract non-unique reads (with poor quality) and duplicated 
@author: irina
'''
import pysam, argparse,sys

def find_nonU(sam, outputs):
    samf_un = pysam.AlignmentFile(outputs[0], "wb", template=sam)
    samf_dup = pysam.AlignmentFile(outputs[1], "wb", template=sam)
    sam_fe = sam.fetch()
    star = 0
    seq = ''
    for read in sam_fe:
        r_star = read.get_reference_positions()[0]
        if r_star == star:
            if read.query_sequence in seq or seq in read.query_sequence:
                samf_dup.write(read)
        else: 
            star = r_star
            seq = read.query_sequence
        if read.mapping_quality <= 3:
            samf_un.write(read)
        

if __name__ == '__main__':
    
    optparser = argparse.ArgumentParser()
    optparser.add_argument('-f', type = str, dest="BamFile", default = "", help = "Bam file")
    optparser.add_argument('-o', type = str, nargs=2, default = ["non_Uniq.bam", "duplicated.bam"], dest="Output", help = "Output filename")
    
    args = optparser.parse_args()
    if len(sys.argv) ==1:
        print optparser.print_help()
        sys.exit(1)
    
    samfile = pysam.AlignmentFile(args.BamFile, "rb")
    find_nonU(samfile, args.Output)