python run_ranked_enrichment.py 
python run_ranked_enrichment.py GSEA --help
python run_enrichment.py  --help

python run_enrichment.py parent-child --help



python run_ranked_enrichment.py  GSEA -o test1.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html  


# cat test1.out

python run_ranked_enrichment.py parent-child -o test2.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html  -m intersection 

# cat test2.out


python run_enrichment.py parent-child -o test3.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html -c bonferroni bh_fdr  -m union 
python run_enrichment.py parent-child -m union -o test4.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html -c bonferroni,bh_fdr 



python run_enrichment.py term-for-term  -o test5.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -c bonferroni -f html  
python run_enrichment.py  term-for-term -o test6.out -i test/deseq.txt -g test/go_test.obo -a test/ga_test.fb  -f html  -c bonferroni