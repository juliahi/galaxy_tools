python run_ranked_enrichment.py -o test1.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html   GSEA 


# cat test1.out

python run_ranked_enrichment.py -o test2.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html parent-child  -m intersection 

# cat test2.out