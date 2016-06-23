#!/usr/bin/python

import sys
import os
import math
import argparse

import Bio.Ontology
import Bio.Ontology.IO as OntoIO



def read_list(filename):
    
    out = []
    with open(filename, 'r') as file_in:
        line = file_in.readline()  #header
        if not (line[0] == '!' or line[0] == '#'): file_in.seek(0)

        for line in file_in:
            content = line.strip().split('\t') 
            if len(content) <= 1:
                print content
                if content[0] != "":
                    out.append(content[0])
            elif content[1] == '1':
                if len(content[0].split('_')) < 2:
                    out.append( content[0])
                else:
                    out.append( "_".join(content[0].split('_')[1:-1]))
            elif content[1] != '0':
                raise Exception("Invalid values in list of genes")
        print out
    return out




def run_term(assocs, go_graph, gene_list, corrections):
    from Bio.Ontology import TermForTermEnrichmentFinder

    ef = TermForTermEnrichmentFinder(assocs, go_graph)
    result = ef.find_enrichment(gene_list, corrections)
    return result

def run_parent_child(assocs, go_graph, gene_list, corrections, method):
    from Bio.Ontology import ParentChildEnrichmentFinder

    ef = ParentChildEnrichmentFinder(assocs, go_graph)
    result = ef.find_enrichment(gene_list, corrections, method)
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
    main_parser = argparse.ArgumentParser(description='run Gene Ontology')
    subparsers = main_parser.add_subparsers(dest='which', help='type of enrichment analysis')
    subparsers.required = True
   
    parser = argparse.ArgumentParser(add_help=False)
    
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-o', '--out', type=str, required=True,
                   help='output file')
    
    required.add_argument('-i', '--inp', type=str, required=True,
                   help='input gene list file')
    required.add_argument('-a', '--assoc', type=str, required=True,
                   help='input associations file (.gaf)')
    required.add_argument('-g', '--gograph', type=str, required=True,
                   help='input GO graph file (.obo)')
    
    
    parser.add_argument('-f', '--outputformat', choices=["html","txt", "gml", "png"],  
                   help='output file format', default = "html")
    
    parser.add_argument('-c', '--corrections', choices=["bonferroni","bh_fdr", "bonferroni,bh_fdr", "bh_fdr,bonferroni"],  
                   help='multiple hypothesis testing corrections', nargs='+', default=[])
    
    
    parser1 = subparsers.add_parser("term-for-term", parents=[parser])
    parser2 = subparsers.add_parser("parent-child", parents=[parser])
    
    

    
    #Term-for-term params
    #parser1.set_defaults(which='term')

    #Parent-child
    parser2.add_argument('-m', '--method', choices=["union", "intersection"],  
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
    
    cors = []
    for cor in args.corrections:
        if "," in cor:
            cors += cor.split(",")
        else:
            cors.append(cor)
    args.corrections = list(set(cors))
    
    
    gene_list = read_list(args.inp)
    
    #gene_rank = [('FBgn0043467', 0.1), ('FBgn0010339', 0.7), ('FBgn0070057', 0.4), ('FBgn0070052', 0.9)]
    
    
    #go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
    #assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")
    
    go_graph = OntoIO.read(args.gograph, "obo")
    assocs = OntoIO.read(args.assoc, "gaf")
    result=None
    
        
    if args.which == "term-for-term":
        result = run_term(assocs, go_graph, gene_list, args.corrections)
    elif args.which == "parent-child":
        #parser.error("Method unimplemented!")
        result = run_parent_child(assocs, go_graph, gene_list, args.corrections, args.method)


    print result
    with open(args.out, 'w+') as outfile:
        assert result!= None,  "An error occured while computing result"
        OntoIO.pretty_print(result, go_graph, outfile, args.outputformat, go_to_url="http://amigo.geneontology.org/amigo/term/")

        
        
        
        
    
    
    
if __name__ == "__main__":
   main()












