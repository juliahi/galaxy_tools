#Filter paired-end fastq files after removing some reads from ends separately

#author: Julia Herman-Izycka
#added: 18.02.2015
#modified: 12.03.2015

import sys, argparse
from strtime import strtime


def writeone(fi, fo, line):
    if fo != None:
        fo.write(line)
        fo.write(fi.readline())
        fo.write(fi.readline())
        fo.write(fi.readline())
    else:
        fi.readline()
        fi.readline()
        fi.readline()
    

def joinPE(m, f1, f2, o1, o2, o3, o4):

    line1=f1.readline()
    line2=f2.readline()
    counter_unsorted = 0
    counter_pairs = 0
    count_left = 0
    count_right = 0
    counter_removed = 0
    while True:
        line=m.readline()
        if line == "": break
        counter_unsorted += 1
        m.readline()
        m.readline()
        m.readline()
        i = line.split()[0] #id of pair
        if line1: i1 = line1.split()[0]
        else: i1 = ''
        if line2: i2 = line2.split()[0]
        else: i2 = ''
        if (i == i1) and (i == i2): 
            writeone(f1, o1, line1)
            line1 = f1.readline()
            writeone(f2, o2, line2)
            line2 = f2.readline()
            counter_pairs+=1
        elif (i == i1):    
            writeone(f1, o3, line1)
            line1 = f1.readline()
            count_left += 1
        elif (i == i2) and (o4 == None):          
            writeone(f2, o3, line2)
            line2 = f2.readline()
            count_right += 1
        elif (i == i2):       
            writeone(f2, o4, line2)
            line2 = f2.readline()
            count_right += 1
        else:
            counter_removed += 1
        #if counter_unsorted % 1000000 == 0:
        #    print "%d paired / %d " % (counter_pairs, counter_unsorted)
    #print """Unfiltered pairs: %d, after filtering: %d pairs,
    #        %d only left, %d only right""" % (counter_unsorted, counter_pairs, count_left, count_right)
    return counter_pairs, count_left+count_right, counter_removed  

def main():
    files = [None]*7
    parser = argparse.ArgumentParser(description='Join paired-end Fastq files after filtering')
    parser.add_argument('-u', type=str, required=True,
                   help='unfiltered fastq')
    parser.add_argument('-i', type=str, nargs=2, required=True,
                   help='filtered fastq left and right')
    
    parser.add_argument('-o', type=str, nargs=2, required=False,
                   help='filtered paired fastq left and right, default:inputname.joined')
    
    group = parser.add_argument_group('singles')
    group.add_argument('-s', type=str, required=False,
                   help='unpaired filtered fastq left and right together')
    group.add_argument('-t',  type=str, nargs=2, required=False,
                   help='unpaired filtered fastq left and right')
    
    args = parser.parse_args()
    print strtime()
    print "Restoring paired-end reads after filtering. "

    files[0] = open(args.u, "r")
    files[1] = open( args.i[0], "r")
    files[2] = open( args.i[1], "r")
    if args.o:
            files[3] = open( args.o[0], "w+")
            files[4] = open( args.o[1], "w+")
    else:
            files[3] = open( args.i[0]+".joined", "w+")
            files[4] = open( args.i[1]+".joined", "w+")
    if args.s:
            files[5] = open( args.s, "w+")
    if args.t:
            files[5] = open( args.t[0], "w+")
            files[6] = open( args.t[1], "w+")
        
     
    number = joinPE(*files)
    
    print "Restored %d pairs, found %d singles, %d pairs removed"% number
    print strtime()
    for f in files:
        if f:
            f.close()
    
    
if __name__ == "__main__":
   main()
