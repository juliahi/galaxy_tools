
import re
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Coverage BedGraph')
    
    parser.add_argument('-i', '--input', type=str, required=True,
                help='input fastq file')
    parser.add_argument('-o', '--output', type=str, required=True,
                help='output file')
    
    parser.add_argument('-p', '--patterns', type=str, nargs='+', required=True,
                help='regular expression to count')
    
    args = parser.parse_args()

    mapped_chars = { '\'' :'__sq__', '\\' : '__backslash__' }


    #with open(args.input, 'r') as inputf:
    fastq_parser = SeqIO.parse(args.input, "fastq") 
    
    n = len(args.patterns)
    counts = [0 for i in xrange(n)]
    
    for key, value in mapped_chars.items():
         for pattern in args.patterns:
               pattern = pattern.replace(value, key)
               
               
    regexps = map(re.compile, args.patterns)
    for record in fastq_parser:
        s = str(record.seq)
        for i in xrange(n):
            if regexps[i].search(s) != None:
                counts[i] += 1

    with open(args.output, 'w') as output:        
        for count, pattern in zip(counts, args.patterns):
            print "Found %d occurrences of pattern %s in file %s"% (count, pattern, args.input)
            output.write("%s\t%d\n"%(pattern, count))

if __name__ == "__main__":
    main()
