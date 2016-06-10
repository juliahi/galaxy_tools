'''
Created on 06-02-2015

@author: irina
modified by julia on 23-02-2015, 21-03-2015
'''
import pysam, time, sys, optparse
import os.path
import os
from strtime import strtime

def is_forward(read, revert):
    if (read.is_read2 and read.is_reverse) or (read.is_read1 and not read.is_reverse): return not revert
    else: return revert



def filtr_bam (f, samf, m_qual, outputs, revert):
    print "All filters, using quality %d"%m_qual
    #samf = pysam.AlignmentFile(f, "rb")
    #samf_for = pysam.AlignmentFile(f.split(".")[0]+"_for.bam", "wb", template=samf)
    #samf_back = pysam.AlignmentFile(f.split(".")[0]+"_back.bam", "wb", template=samf)
    samf_for = pysam.AlignmentFile(outputs[0], "wb", template=samf)
    samf_back = pysam.AlignmentFile(outputs[1], "wb", template=samf)
    sam_fe = samf.fetch()
    counter = 0
    size = 0
    for read in sam_fe:
        size += 1
        if read.is_secondary or read.is_unmapped or read.mapping_quality <= m_qual: continue
        #=======================================================================
        # elif read.is_unmapped: 
        #     continue
        # elif read.mapping_quality <= 3: continue
        #=======================================================================
        elif is_forward(read, revert):
        #elif read.is_reverse:
            samf_for.write(read)
        else: samf_back.write(read)
        counter += 1
    print "%d out of %d read pairs remained"%(counter, size)
    samf_for.close()
    samf_back.close()
        
def strand_sep(f, samf, outputs, revert):
    print "Strand separation"
    #samf_for = pysam.AlignmentFile(f.split(".")[0]+"_for.bam", "wb", template=samf)
    #samf_back = pysam.AlignmentFile(f.split(".")[0]+"_back.bam", "wb", template=samf)
    samf_for = pysam.AlignmentFile(outputs[0], "wb", template=samf)
    samf_back = pysam.AlignmentFile(outputs[1], "wb", template=samf)
    sam_fe = samf.fetch()
    counter = 0
    for read in sam_fe:
        if is_forward(read, revert):
            samf_for.write(read)
        else: samf_back.write(read)
        counter += 1
    print "%d read pairs remained"%counter 
    
    samf_for.close()
    samf_back.close()

def only_qual(f, samf, outputs, revert, m_qual, ag=False, m=False, q=False, al=False, sep=False):
    print "Quality filters"
    sam_fe = samf.fetch()
    counter = 0
    size = 0
    if sep:
        print "with strand separation"
        #samf_for = pysam.AlignmentFile(f.split(".")[0]+"_for.bam", "wb", template=samf)
        #samf_back = pysam.AlignmentFile(f.split(".")[0]+"_back.bam", "wb", template=samf)
        samf_for = pysam.AlignmentFile(outputs[0], "wb", template=samf)
        samf_back = pysam.AlignmentFile(outputs[1], "wb", template=samf)
    else:
        #sam_clean = pysam.AlignmentFile(f.split(".")[0]+"_clean.bam", "wb", template=samf)
        sam_clean = pysam.AlignmentFile(outputs[0], "wb", template=samf)
    if ag and not m and not q:
        print "Filtering first_alignment"
        for read in sam_fe:
            size += 1
            if read.is_secondary: continue 
            if sep:
                if is_forward(read, revert):
                    samf_for.write(read)
                else: samf_back.write(read)
            else: sam_clean.write(read)
            counter += 1
    elif m and not ag and not q:
        print "Filtering mapped_reads"
        for read in sam_fe:
            size += 1
            if read.is_unmapped: continue
            if sep:
                if is_forward(read, revert):
                    samf_for.write(read)
                else: samf_back.write(read)
            else: sam_clean.write(read)
            counter += 1
    elif q and not m and not ag:
        print "Filtering higher_quality, using quality %d"%m_qual
        for read in sam_fe:
            size += 1
            if read.mapping_quality <= m_qual: continue
            if sep:
                if is_forward(read, revert):
                    samf_for.write(read)
                else: samf_back.write(read)
            else: sam_clean.write(read)
            counter += 1
    elif q and m and not ag:
        print "Filter higher_quality and mapped_reads, using quality %d"%m_qual
        for read in sam_fe:
            size += 1
            if read.mapping_quality <= m_qual or read.is_unmapped: continue
            if sep:
                if is_forward(read, revert):
                    samf_for.write(read)
                else: samf_back.write(read)
            else: sam_clean.write(read)
            counter += 1
    elif q and ag and not m:
        print "Filter higher_quality and alignment, using quality %d" %m_qual
        for read in sam_fe:
            size += 1
            if read.mapping_quality <= m_qual or read.is_secondary: continue
            if sep:
                if is_forward(read, revert):
                    samf_for.write(read)
                else: samf_back.write(read)
            else: sam_clean.write(read)
            counter += 1
    elif m and ag and not q:
        print "Filter mapped_reads and alignment"
        for read in sam_fe:
            size += 1
            if read.is_unmapped or read.is_secondary: continue
            if sep:
                if is_forward(read, revert):
                    samf_for.write(read)
                else: samf_back.write(read)
            else: sam_clean.write(read)
            counter += 1
    elif al:
        print "All filters, using quality %d"% m_qual
        for read in sam_fe:
            size += 1
            if read.is_secondary or read.is_unmapped or read.mapping_quality <= m_qual: continue 
            if sep:
                if is_forward(read, revert):
                    samf_for.write(read)
                else: samf_back.write(read)
            else: sam_clean.write(read)
            counter += 1
    else: print "Option is unknown", ag, m, q, al
    print "%d out of %d read pairs remained"%(counter, size) 
    if sep:
        samf_for.close()
        samf_back.close()
    else:
        sam_clean.close()
        
        
def paring(sam):
    sam_fe = sam.fetch()
    total_time = time.time() + 3
    names = {}
    r_num = 0
    r_bad = 0
    for r in sam_fe:
        r_num +=1
        if r.qname in names.keys(): 
            print r, "\n",  names[r.qname][:]
            r_bad +=1
            names[r.qname].append(r)
        else: names[r.qname]=[r]
        try: 
            pysam.AlignmentFile.mate(sam,r)
            print  r.qname, "kolega", pysam.AlignmentFile.mate(sam, r).qname
        except ValueError as e: pass
    zle = float(r_bad)/float(r_num) *100
    print zle
            #print "Problem", str(e).split()[-2:] 
        #if time.time() > total_time: sys.exit(1)
    
    
if __name__ == '__main__':
    
    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-f', type = "string", default = "", dest="BamFile", help = "Bam file")
    optparser.add_option('-i', type = "string", default = "", dest="IndexFile", help = "Index file")
    optparser.add_option('-o', type = "string", nargs=2, default = ["for.bam", "back.bam"], dest="Output", help = "Output filename")
    optparser.add_option('-a', dest="AllFilters", action = "store_true", default = False, help = "All filters will used: only first alignment, only mapped reads, with mapping quality > 3 or given value with -m option, each strand in separate file")
    
    optparser.add_option('-q', dest="OnlyQuality", action = "store_true", default = False, help = " All quality filters will be run: first_afignment, mapped_reads, reads with higher, than given with option -m quality")
    optparser.add_option('-g', dest="Alignment", action = "store_true", default = False, help = "Running first_alignment filter")
    optparser.add_option('-p', dest="Mapped", action = "store_true", default = False, help = "Filter reads that are mapped")
    optparser.add_option('-l', dest="MQuality", action = "store_true", default = False, help = "Filter reads with higher than given with option -m mapping quality")
    
    optparser.add_option('-m', type = "int", default = 3, dest="MapQ", help = "Mapping quality that will be used with -a and -qa options. Default = 3")
    optparser.add_option('-s', dest="StrandSep", action = "store_true", default = False, help = "Strands into two files")
    optparser.add_option('-r', dest="Revert", action = "store_true", default = False, help = "Use second read as strand indicator")
    
    (opts,args) = optparser.parse_args()
    if len(sys.argv) ==1:
        print optparser.format_help()
        sys.exit(1)
        
    if opts.AllFilters and opts.OnlyQuality:
        optparser.error("Options -a and -q are mutually exclusive")
        
    if not (opts.AllFilters or opts.OnlyQuality or opts.Alignment or opts.Mapped or opts.MQuality or opts.StrandSep):
        optparser.error("Chose one of following options: -a, -q, -g, -p, -l, -s")
    
    if (opts.OnlyQuality  and opts.Alignment) or (opts.OnlyQuality  and opts.Mapped) or (opts.OnlyQuality  and opts.MQuality ):   
        optparser.error("Option -q and -g as well as -q and -p or -q and -l are mutually exclusive")
        
    #if opts.Alignment and opts.Mapped and opts.MQuality:   
    #    optparser.error("Run with option -q instead")    
        
    #if (opts.AllFilters and opts.StrandSep) or (opts.OnlyQuality and  opts.StrandSep):
    #    optparser.error("Options -a and -s as well as -q and -s are mutually exclusive")
    
    outputs=opts.Output
    print strtime()

    bam_file = opts.BamFile 
    if opts.IndexFile=="":
        pysam.index(bam_file)
        print strtime() + ": Index created"
    elif not os.path.isfile(bam_file+".bai"):
        try:
            os.symlink(opts.IndexFile, bam_file+".bai")
        except Exception:
            print "Cannot link index file %s. Check if file exists" % opts.IndexFile
            
            exit(1)
        print "Index linked"

    samfile = pysam.AlignmentFile(bam_file, "rb")
        
    #samsamfile_fe = samfile.fetch()
    if opts.AllFilters:
        filtr_bam(bam_file, samfile, opts.MapQ, outputs, opts.Revert)
    elif opts.Alignment and opts.Mapped and opts.MQuality and opts.StrandSep:
        filtr_bam(bam_file, samfile, opts.MapQ, outputs, opts.Revert)
    elif opts.Alignment and opts.Mapped and opts.MQuality:
        only_qual(bam_file, samfile, outputs, opts.Revert, opts.MapQ, False, False, False, True, sep=opts.StrandSep)
    elif opts.OnlyQuality or opts.Alignment or opts.Mapped or opts.MQuality :
        only_qual(bam_file, samfile, outputs, opts.Revert, opts.MapQ, opts.Alignment, opts.Mapped, opts.MQuality, opts.OnlyQuality, sep=opts.StrandSep)
    elif opts.StrandSep:
        strand_sep(bam_file, samfile, outputs, opts.Revert)
        
    else: print ""
    print strtime()+": Finished"
        
    #paring(samfile)
    
    
