# 2. Mapping reads to intergenic regions 

First, we can use the following python script to generate intergenic GFF/GTF files, providing the _B. coprophila_ genome annotation (-g), along with metrics to describe the minimum distance from genes and the minimum/maximum length for the intergenic regions. For this study we ran the script like so: 
```
python3 01_get_intergenic_GFF3.py
    -g ../../Annotations/bcop_core_GRC.gff3
    --min_distance 500 --min_length 1000
    --max_length 20000
```
Output (intergenic GFF and GTF) from 01_get_intergenic_GFF3.py can be fouind in /outputs/ 

Then we can use `02_intergenic_TPM.sh`, which runs StringTie on the uniquley mapped BAM files (generated in step 01), this time providing the intergenic gtf to calculate intergenic TPMs
```
for file in $(ls *_uniquely_mapped.bam)
do
    base=$(basename "$file" "_uniquely_mapped.bam")
    output_gtf="${base}_intergenic.gtf"
    echo "Calculating TPM for $file"
    stringtie "$file" -p 16 -o "$output_gtf" -G bcop_core_GRC.intergenic.gtf
done
```
Then we can run another custom python script `03_get_intergenic_TPM.py` in the same directory where we ran StringTie to generate a combined TPM file 
```
python3 03_get_intergenic_TPM.py -t . 
```
Finally, we can run the R script `04_TPM_cutoff_intergenic_deconvolution.R`. This script classifies genes as expressed on not expressed. See paper methods and Supplementary materials for detailed method

> **Note** This script requires both 
> `Bcop_GRC_quiescence/01_RNAseq_mapping/outputs/combined_TPM_only.tsv` as well as
> `Bcop_GRC_quiescence/02_intergenic_TPM/outputs/combined_intergenic_TPM.tsv`

This visualises and performs statistics to generate a TPM cutoff value for expression, based on this intergenic mapping rate. The final output is _TPM_genes.tsv_ 

| Output file | Description |
|-------------|-------------|
| `TPM_genes.tsv` | Final output classifying all _B. coprophila_ genes as either Expressed on Not Expressed based |
| `gene_overlap_0-4h_TPM.tsv` | All GRC-linked genes with TPM > 0 in 0-4h embryos libraries |
| `gene_overlap_4-8h_TPM.tsv` | All GRC-linked genes with TPM > 0 in 4-8h embryos libraries |
| `gene_overlap_late-larva-early-pupa_TPM.tsv` | All GRC-linked genes with TPM > 0 in late-larva/early pupae libraries |
| `gene_overlap_adult_TPM.tsv` | All GRC-linked genes with TPM > 0 in adult libraries |
