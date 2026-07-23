# 3. Defining confidently expressed GRC-linked genes and downstream analysis 

For a GRC-linked gene to be considered expressed in this study, it must satisfy a few conditions:
- Have an TPM > cutoff generated in /02_intergenic_TPM/ in at least 2/3rds of replicate libraries (per sex per development stage)
- See no significant expression in somatic libraries 
- Not be highly similar core-chromosome paralogue (i.e. no BLASTn paralogue with a % identity * % coverage > 70) 
 
Hence these scripts mainly sort the data to identify confidently expressed GRC-linked genes in this manner, before performing downstream, analysis 

First, loading the _TPM_genes.tsv_ (Bcop_GRC_quiescence/02_intergenic_TPM/outputs/TPM_genes.tsv) into R, we run:
```
01_classify_GRC_expressed_genes.R 
```
This give us the following useful output tables:
1. _Expressed_GRC_genes.tsv_ (full table of all genes considered expressed)
2. _Expressed_GRC_genes_summary.tsv_ (summarised table)
3. _Corrected_expressed_GRC_genes.tsv_ table (whereby false positive genes with significant smaitc mismapping have been removed)

Next we can then use the _Corrected_expressed_GRC_genes.tsv_ output to retrieve nucelotide sequences for confidently expressed GRC-linked genes like so
```
python 02_get_GRC_core_fastas.py -g ../Annotations/idBraCopr2.1.primary.masked a ../Annotations/bcop_core_GRC.gff3 -s Corrected_expressed_GRC_genes.tsv
```

