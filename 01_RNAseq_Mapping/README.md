# 1. RNAseq mapping code

This script maps RNA reads from embryo, larval and adult data and maps it to the refernece _B. coprophila_ gneome and annotation 

First reads are trimmed (QC had been done prior to this study)
```
echo "Trimming reads with fastp..."
for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
	echo "Trimming $base"
  	fastp -i ${base}_1.fq.gz -I ${base}_2.fq.gz -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz && \
  	rm ${base}_1.fq.gz ${base}_2.fq.gz
done
```
Next we index the _B. coprophila gneome_ using STAR 
```
echo "run genomeGenerate"
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeSAindexNbases 12 \
--outFileNamePrefix idBraCopr2.1.primary.masked.fa \
--sjdbGTFfile idBraCopr2.1.primary.masked_core_and_grc_braker3.gff3  \
--genomeDir idBraCopr2.1.primary.masked.fa.STAR \
--genomeFastaFiles idBraCopr2.1.primary.masked.fa
```
Then we map the reads to the genome using STAR, with filters to ensure proper read mapping and providing the annotation to increase fidelity 
```
for file in $(ls *_1.trimmed.fq.gz)
do
    base=$(basename "$file" "_1.trimmed.fq.gz")
    echo "Aligning $base to B_coprophila genome"
    STAR \
        --runThreadN 16 \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --readFilesIn ${base}_1.trimmed.fq.gz ${base}_2.trimmed.fq.gz \
        --outTmpDir ${base}.out \
        --outFileNamePrefix ${base}.STAR. \
        --outFilterMismatchNmax 2 \
        --genomeDir idBraCopr2.1.primary.masked.fa.STAR
done
```
Then we apply post-mapping filtering for uniquely mapped reads (MAPQ=255) and run StringTie
```
for file in $(ls *.bam)
do
    base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
    output_file=${base}_uniquely_mapped.bam
    echo "Filtering unique reads for $file"
    samtools view -@ 16 -q 255 -b $file > $output_file
done

samtools index -M *_uniquely_mapped.bam

for file in $(ls *_uniquely_mapped.bam)
do
    base=$(basename $file "_uniquely_mapped.bam")
    output_file=${base}_unique_mapping_TPM.gtf
    echo "Calculating TPM for $file"
    stringtie -o $output_file -G idBraCopr2.1.primary.masked_core_and_grc_braker3.gff3  $file
done
```
Finally, we can pass the StringTie outputs directly to our custom script simply by running
```
python3 02_get_TPM_values.py -t .
```
Outputs from 02_get_TPM_values.py can be found in /outputs/
