
#!/usr/bin/python3


import sys;
import re;


def use_regexp(element, s):
    expr = element+"([\s=:]\"?[\w\.]+\"?)";
    transid=re.findall(expr, s);
    
    if len(transid) == 0: 
        return ''
    return transid[-1].strip("=\" :") 

def printbedline(estart,eend,field,nline):
    try:
        # use regular expression to get transcript_id, gene_id and expression level
        geneid=use_regexp('gene_id',field[8])
        genename=use_regexp('gene_name',field[8])
        transid=use_regexp('transcript_id',field[8])
        
        #if len(geneid)==0:
        #    print('Warning: no gene_id field ',file=sys.stderr);
        #if len(transid)==0:
        #    print('Warning: no transcript_id field',file=sys.stderr);
        #    transid='Trans_'+str(nline);
        if transid in allids.keys():
            transid2=transid+'_DUP'+str(allids[transid]);
            allids[transid]=allids[transid]+1;
            transid=transid2;
        else:
            allids[transid]=1;

        
        perm = sorted(range(len(estart)), key=lambda k: estart[k])
        estart = [str(estart[k]) for k in perm]
        eend = [str(eend[k]) for k in perm]

        starts=','.join(estart);
        ends=','.join(eend);
        
        chrom = field[0]
        #if not chrom.startswith('chr'):
        #    if chrom=='MT': chrom='M'
        #    elif len(chrom) > 2: chrom='_'.join(chrom.lower().split('.')[::-1])
        #    chrom = 'chr'+chrom
        strand = field[6]
       
        if genename:
            name = genename
        else:
            name = use_regexp('gene', field[8])
        if geneid == '': geneid =use_regexp('GeneID', field[8]) 
        output.write(transid+'\t' + chrom + "\t" + strand + "\t" + str(estart[0]) + "\t" + str(eend[-1]) + "\t" + starts+'\t'+ends + "\t" + name + "\t" + geneid + "\n" );

    except ValueError:
        print('Error: non-number fields at line '+str(nline),file=sys.stderr);





if __name__ == "__main__":

    if len(sys.argv)!=3:
        print('Usage: [GTF file] [output_file] \n');
        print('\nNote:');
        print('1\tOnly "exon" and "transcript" are recognized in the feature field (3rd field).');
        print('2\tIn the attribute list of .GTF file, the script tries to find "gene_id", "transcript_id" and "FPKM" attribute, and convert them as name and score field in .BED file.');
        print('Author: Wei Li (li.david.wei AT gmail.com), modified by Julia H.-I.');
        sys.exit();


    allids={};
    estart=[];
    eend=[];

    nline=0;

    prevfield=[];
    prevtransid=None;
    output = open(sys.argv[2], 'w+')

    for lines in open(sys.argv[1]):
        if lines[0] == '#': continue
        field=lines.strip().split('\t');
        nline=nline+1;

        if len(field)<9:
            print('Error: the GTF should has at least 9 fields at line '+str(nline),file=sys.stderr);
            continue;

        if field[2]!='exon' and field[2] !='transcript':
            continue;
        #remove unlocalized genomic contigs
        if len(field[0]) > 5 and field[0].startswith('chr'):
            continue;

        transid=use_regexp('transcript_id',field[8]);
        if field[2]=='transcript' or (transid!='' and transid != prevtransid):
                #print('prev:'+prevtransid+', current:'+transid);
                # A new transcript record, write
                if len(estart)!=0:
                    printbedline(estart,eend,prevfield,nline);
                    
                estart=[];
                eend=[];
                prevfield=field;
                prevtransid=transid;

        if 'exon' in field[2]:
            try:
                        est=int(field[3])-1;
                        eed=int(field[4]);
                        if est > eed:
                            est, eed = eed, est

                        estart+=[est];
                        eend+=[eed];
            except ValueError:
                        print('Error: non-number fields at line '+str(nline),file=sys.stderr);


    # the last record
    if len(estart)!=0 and len(prevfield)>8:
        printbedline(estart,eend,prevfield,nline);
        
    output.close()

