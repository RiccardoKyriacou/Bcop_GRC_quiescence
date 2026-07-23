# 3. Defining confidently expressed GRC-linked genes and downstream analysis 

First, loading the _TPM_genes.tsv_ (Bcop_GRC_quiescence/02_intergenic_TPM/outputs/TPM_genes.tsv) into R, we run:
```
01_classify_GRC_expressed_genes.R 
```
This give us useful output tables _Expressed_GRC_genes.tsv_ (full table of all genes considered expressed), _Expressed_GRC_genes_summary.tsv_ (summarised table) and the final _Corrected_expressed_GRC_genes.tsv_ table (whereby false positive genes with significant smaitc mismapping have been removed)
