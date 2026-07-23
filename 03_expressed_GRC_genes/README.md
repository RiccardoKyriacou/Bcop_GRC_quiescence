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
This give us the following useful output tables (03_expressed_GRC_genes/outputs/expressed_GRC_gene_tables):
1. _Expressed_GRC_genes.tsv_ (full table of all genes considered expressed)
2. _Expressed_GRC_genes_summary.tsv_ (summarised table)
3. _Corrected_expressed_GRC_genes.tsv_ table (whereby false positive genes with significant smaitc mismapping have been removed)

Next we can then use the _Corrected_expressed_GRC_genes.tsv_ output to retrieve nucelotide sequences for both the confidently expressed GRC-linked genes, as well as all core chromosome genes in _B.coprophila_. We run the script like so
```
python 02_get_GRC_core_fastas.py -g ../Annotations/idBraCopr2.1.primary.masked a ../Annotations/bcop_core_GRC.gff3 -s Corrected_expressed_GRC_genes.tsv
```

We can then run the following bash script
```
03_BLAST_GRCgenes.sh
```
This allows us to perform a BLASTn search of the expressed GRC-linked gene sequences we just generated. We BLASTn these queries first against the complete list of annotated core chromosome genes (generated in the previous step), and then the full core reference genome. This is done to identify GRC-gene similarity to core-chromosome paralogues.
```
### STEP 1) BLAST GRC genes against core genes ###

# Make a blast data base of core-chromosome genes 
echo "making BLASTDB from core chromosome genes"
makeblastdb \
-in core_genes.fasta \
-dbtype nucl \
-parse_seqids \
-out core_genes_DB

# Blast GRC genes against core genes 
echo "BLASTn for GRC genes"
blastn \
-query GRC_genes_to_BLAST.fasta \
-db core_genes_DB \
-out GRC_gene_BLAST_output.tsv \
-outfmt '6 std qlen slen qseq sseq'
```
When we BLAST the GRC-linked gene against the core genome (in case GRC-gene paralogues are not annotated correctly), we utilise a python script to mask out the GRC genome (located in /outputs/) 
```
### STEP 2) BLAST GRC genes against masked core genome ###

# Sync files in
# GRC genome
rsync -av /mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked.fa $SCRATCH
# Masking script
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST/scripts/mask_genome.py $SCRATCH

# Make masked GRC genome
python3 mask_genome.py -g idBraCopr2.1.primary.masked.fa -o masked_GRC_genome.fasta -m GRC

# Make a blast data base of core-chromosome genes 
echo "making BLASTDB from maksed genome"
makeblastdb \
-in masked_GRC_genome.fasta \
-dbtype nucl \
-parse_seqids \
-out masked_GRC_genome_DB

# Blast GRC genes against whole genome 
echo "BLASTn for GRC genes"
blastn \
-query GRC_genes_to_BLAST.fasta \
-db masked_GRC_genome_DB \
-out GRC_genome_BLAST_output.tsv \
-outfmt '6 std qlen slen qseq sseq'
```
These BLAST results (located in outputs/GRC_v_Core_BLAST/) are then passed through another custom python script, 03_get_BLAST_table.py, in order to generate a table of expressed GRC genes and their respective best-hit core paralogue. This allows us to identify any GRC-linked genes with highly similar core paralogues which decreasue our condifence in their expression.

We pass _Corrected_expressed_GRC_genes.tsv_ generated in step 01 (-t), the GRC vs core gene (-c) and the GRC vs core genome outputs (-g) like so 
```
python 04_get_BLAST_table.py -t GRC_gene_expression.tsv -c GRC_v_Core_gene_BLAST_output - g GRC_v_Core_genome_BLAST_output.tsv
```
Our analysis now moves onto focusing on the confidently expressed GRC genes. 

