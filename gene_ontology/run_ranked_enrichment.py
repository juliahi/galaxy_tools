#!/usr/bin/python

import sys
import os
import math
import argparse

import Bio.Ontology
import Bio.Ontology.IO as OntoIO



def read_deseq_output(filename):
    column=5 #value from which column to take -> 5=log2FoldChange
    remove_inf = True   #remove entries with infinity? (otherwise change to max)
    
    out = []
    with open(filename, 'r') as file_in:
        file_in.readline()  #header

        maxval = 0
        minval = 0

        for line in file_in:
            content = line.split('\t') 
            if len(content[0].split('_')) < 2:
                print "problem in line:", line
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
    parser = argparse.ArgumentParser(description='run Ranked Gene Ontology on DESeq output')
    
    
    parser.add_argument('-o', '--out', type=str, required=True,
                   help='output file')
    
    parser.add_argument('-i', '--inp', type=str, required=True,
                   help='input DESeq result file')
    parser.add_argument('-a', '--assoc', type=str, required=True,
                   help='input associations file (.gaf)')
    parser.add_argument('-g', '--gograph', type=str, required=True,
                   help='input GO graph file (.obo)')
    
    
    parser.add_argument('-t', '--type', choices=["GSEA", "parent-child"],
                   help='type of enrichment analysis', required=True)
    
    parser.add_argument('-f', '--outputformat', choices=["html","txt", "gml", "png"],  
                   help='output file format', default = "html")
    
    
    #GSEA params
    parser.add_argument('-p', '--perms', type=int, default=1000, 
                   help='number of permutations to compute p-values')
    parser.add_argument('-n', '--minset', type=int, default=1, 
                   help='minimal intersection between set of genes in rank and in studied set')
    parser.add_argument('-c', '--corrpower', type=float, default=1, 
                   help='how strong correlation will affect enrichment')


    
    #validate args
    args = parser.parse_args()
    check_file(parser, args.inp, 'r')
    check_file(parser, args.assoc, 'r')
    check_file(parser, args.gograph, 'r')
    check_file(parser, args.out, 'w+')
    
    
    
    if args.type == "GSEA":
        #check parameters
        #TODO
        pass;
    
    gene_rank = read_deseq_output(args.inp)
    
    #gene_rank = [('FBgn0043467', 0.1), ('FBgn0010339', 0.7), ('FBgn0070057', 0.4), ('FBgn0070052', 0.9)]
    #go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
    #assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")
    
    print gene_rank
    
    go_graph = OntoIO.read(args.gograph, "obo")
    assocs = OntoIO.read(args.assoc, "gaf")
    result=None
    
    if args.type == "GSEA":
        result = run_gsea(assocs, go_graph, gene_rank, args.perms, args.minset, args.corrpower)
    else:
        parser.error("Method unimplemented!")
        #result = run_parent_child(assocs, go_graph, gene_rank, )
    print result
    

    with open(args.out, 'w+') as outfile:
        assert result!= None,  "An error occured while computing result"
        OntoIO.pretty_print(result, go_graph, outfile, args.outputformat)

        
        
        
        
    
    
    
if __name__ == "__main__":
   main()












