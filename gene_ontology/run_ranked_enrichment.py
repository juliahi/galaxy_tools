#!/usr/bin/python

import sys
import os
import math
import argparse

import Bio.Ontology
import Bio.Ontology.IO as OntoIO


def read_deseq_output(filename, column):
    ##column=5 #value from which column to take -> 5=log2FoldChange
    remove_inf = True   #remove entries with infinity? (otherwise change to max)
    
    out = []
    with open(filename, 'r') as file_in:
        file_in.readline()  #header

        maxval = 0
        minval = 0

        for line in file_in:
            content = line.split('\t') 
            if len(content) <= column:
                continue
            if len(content[0].split('_')) < 2:
                out.append( (content[0], content[column]))
            else:
                out.append( ("_".join(content[0].split('_')[1:-1]), content[column]))
            if not remove_inf and not "Inf" in content[column]:
                maxval = max(maxval, float(content[column]))
                minval = min(minval, float(content[column]))
            
    if remove_inf:
            return [(x[0], float(x[1])) for x in out if "Inf" not in x[1] ]
    else:
            for i, x in enumerate(out):
                if x[1] == "Inf":
                    out[i][1] = maxval
                elif x[1] == "-Inf":
                    out[i][1] = minvalue
                else:
                    out[i][1] = float(x[1])
    return out






def run_gsea(assocs, go_graph, gene_rank, perms, minset, corr):
    from Bio.Ontology import GseaEnrichmentFinder

    ef = GseaEnrichmentFinder(assocs, go_graph)
    result = ef.find_enrichment(gene_rank, perms_no = perms, min_set_rank_intersection=minset, corr_power=corr)

    return result

def run_parent_child(assocs, go_graph, gene_rank, side, corrections, rank_as_population, method):
    from Bio.Ontology import RankedParentChildEnrichmentFinder

    ef = RankedParentChildEnrichmentFinder(assocs, go_graph)
    result = ef.find_enrichment(gene_rank, side, corrections, rank_as_population, method)
    return result




def check_file(parser, arg, openparam):
    if openparam == 'r':
        if not os.path.exists(arg):
            parser.error("The file %s does not exist!" % arg)
    else:
        try:
            f=open(arg, openparam)
            f.close()
        except:
            parser.error("Cannot create file %s" % arg)




def main():
    main_parser = argparse.ArgumentParser(description='run Ranked Gene Ontology on DESeq output')
    subparsers = main_parser.add_subparsers(dest='which', help='type of enrichment analysis')
    subparsers.required = True
   
    parser = argparse.ArgumentParser(add_help=False)
    
    
    
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-o', '--out', type=str, required=True,
                   help='output file')
    
    required.add_argument('-i', '--inp', type=str, required=True,
                   help='input DESeq result file')
    required.add_argument('-a', '--assoc', type=str, required=True,
                   help='input associations file (.gaf)')
    required.add_argument('-g', '--gograph', type=str, required=True,
                   help='input GO graph file (.obo)')
    
    
    parser.add_argument('-f', '--outputformat', choices=["html","txt", "gml", "png"],  
                   help='output file format', default = "html")
    
    parser.add_argument('-x', '--field', type=int, default=6,
                   help='Field from file for ranking genes (starting from 1)')
    
    
    parser1 = subparsers.add_parser("GSEA", parents=[parser])
    parser2 = subparsers.add_parser("parent-child", parents=[parser])
    
    #GSEA params
    parser1.add_argument('-p', '--perms', type=int, default=1000, 
                   help='number of permutations to compute p-values')
    parser1.add_argument('-n', '--minset', type=int, default=1, 
                   help='minimal intersection between set of genes in rank and in studied set')
    parser1.add_argument('-c', '--corrpower', type=float, default=1, 
                   help='how strong correlation will affect enrichment')

    #parser1.set_defaults(which='GSEA')

    #Ranked Parent-child
    parser2.add_argument('-s', '--side', choices=["+","-", "+/-"],  
                   help='from which side (highest or lowest) to start computation', default = "+")
    parser2.add_argument('-c', '--corrections', choices=["bonferroni","bh_fdr"],  
                   help='multiple hypothesis testing corrections', nargs='+', default=[])
    parser2.add_argument('-r', '--rank-as-population', action='store_true',
                   help='only the rank should be used as population')
    parser2.add_argument('-m', '--method', choices=["union", "intersection", "term"],  
                   help='method used to compute probabilities', default = "union")
    #parser2.set_defaults(which='parent-child')
    
    #validate args    
    if len(sys.argv) < 2:
        main_parser.print_usage()
        sys.exit(1)
        
    args = main_parser.parse_args()
    check_file(main_parser, args.inp, 'r')
    check_file(main_parser, args.assoc, 'r')
    check_file(main_parser, args.gograph, 'r')
    check_file(main_parser, args.out, 'w+')
    
    
    
    if args.which == "GSEA":
        #check parameters
        
        if args.perms < 1:
            parser.error('wrong number of permutations: %d' % args.p)
        
    #elif args.which == "parent-child":
    if not 2 <= args.field <= 8 :
        parser.error("Field must be a number from 2 to 8")
    

    gene_rank = read_deseq_output(args.inp, args.field-1)
    
    #gene_rank = [('FBgn0043467', 0.1), ('FBgn0010339', 0.7), ('FBgn0070057', 0.4), ('FBgn0070052', 0.9)]
    
    
    #go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
    #assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")
    
    go_graph = OntoIO.read(args.gograph, "obo")
    assocs = OntoIO.read(args.assoc, "gaf")
    result=None
    
    if args.which == "GSEA":
        result = run_gsea(assocs, go_graph, gene_rank, args.perms, args.minset, args.corrpower)
    else:
        #parser.error("Method unimplemented!")
        result = run_parent_child(assocs, go_graph, gene_rank, args.side, args.corrections, args.rank_as_population, args.method)


    print result
    with open(args.out, 'w+') as outfile:
        assert result!= None,  "An error occured while computing result"
        OntoIO.pretty_print(result, go_graph, outfile, args.outputformat, go_to_url="http://amigo.geneontology.org/amigo/term/")

        
        
        
        
    
    
    
if __name__ == "__main__":
   main()












