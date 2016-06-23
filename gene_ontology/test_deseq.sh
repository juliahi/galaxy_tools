python run_ranked_enrichment_deseq.py 
python run_ranked_enrichment_deseq.py GSEA --help
python run_enrichment_deseq.py  --help

python run_enrichment_deseq.py parent-child --help



python run_ranked_enrichment_deseq.py  GSEA -o test1.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html  


# cat test1.out

python run_ranked_enrichment_deseq.py parent-child -o test2.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html  -m intersection 

# cat test2.out


python run_enrichment_deseq.py parent-child -o test3.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html -c bonferroni bh_fdr  -m union 
python run_enrichment_deseq.py parent-child -m union -o test4.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html -c bonferroni,bh_fdr 



python run_enrichment_deseq.py term-for-term  -o test5.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -c bonferroni -f html  
python run_enrichment_deseq.py  term-for-term -o test6.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html  -c bonferroni
