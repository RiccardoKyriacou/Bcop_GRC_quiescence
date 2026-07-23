# 2. Mapping reads to intergenic regions 

First, we can use the following python script to generate intergenic GFF/GTF files, providing the _B. coprophila_ genome annotation (-g), along with metrics to describe the minimum distance from genes and the minimum/maximum length for the intergenic regions. For this study we ran the script like so: 
```
python3 01_get_intergenic_GFF3.py -g ../../Annotations/bcop_core_GRC.gff3 --min_distance 500 --min_length 1000 --max_length 20000
```
Output (intergenic GFF and GTF) from 01_get_intergenic_GFF3.py can be fouind in /outputs/ 

Then we can use 02_intergenic_TPM.sh, which runs StringTie on the uniquley mapped BAM files (generated in step 01), this time providing the intergenic gtf to calculate intergenic TPMs
```
for file in $(ls *_uniquely_mapped.bam)
do
    base=$(basename "$file" "_uniquely_mapped.bam")
    output_gtf="${base}_intergenic.gtf"
    echo "Calculating TPM for $file"
    stringtie "$file" -p 16 -o "$output_gtf" -G bcop_core_GRC.intergenic.gtf
done
```
Then we can run another custom python script in the same directory where we ran StringTie to generate a combined TPM file 
```
pyhton3 03_get_intergenic_TPM.py -t . 
```
Finally, the R script 04_intergenic_TPM_deconvolution.R visualises and performs statistics to generate a TPM cutoff value for expression, based on this intergenic mapping rate
