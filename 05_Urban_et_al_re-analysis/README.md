# 5 Re-analysis of Urban et al. (2021) pooled embryo data

We first download the pooled embryo data using the script 01_download_pooled_embryo.sh

```
srx_list=(SRX6716708 SRX6716709 SRX6716710 SRX6716711 SRX6716704 SRX6716705)

outdir="embryo_fastqs"
mkdir -p "$outdir"

for srx in "${srx_list[@]}"; do
  echo "Fetching runs for $srx ..."
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srx}&result=read_run&fields=run_accession,fastq_ftp" \
    | tail -n +2 | while IFS=$'\t' read -r run fastqs; do
        IFS=';' read -ra files <<< "$fastqs"
        for url in "${files[@]}"; do
          fname="$outdir/$(basename "$url")"
          if [[ -s "$fname" ]]; then
            echo "  - $fname exists, skipping."
          else
            echo "  - downloading $url"
            wget -c -O "$fname" "ftp://$url"
          fi
        done
      done
done
```
We then use the script 02_STAR_Stringtie.sh to first trim the reads using fastp
```
do
	base=$(basename $file "_1.fq.gz")
	echo "Trimming $base"
  	fastp -i ${base}_1.fq.gz -I ${base}_2.fq.gz -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz && \
  	rm ${base}_1.fq.gz ${base}_2.fq.gz
done
```
Then we use STAR to index the reference _B. coprophila_ genome (providing annotation) and map the RNA-reads. We use the same mapping perameters as in /01_RNAseq_Mapping/01_STAR_TPM_Nmax2.
```
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeSAindexNbases 12 \
--outFileNamePrefix idBraCopr2.1.primary.masked.fa \
--sjdbGTFfile idBraCopr2.1.primary.masked_core_and_grc_braker3.gff3  \
--genomeDir idBraCopr2.1.primary.masked.fa.STAR \
--genomeFastaFiles idBraCopr2.1.primary.masked.fa

echo "Aligning B_coprophila RNAseq reads..."
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
We then filter for only uniquely mapping reads using samtools and run StringTie
```
for file in $(ls *.bam)
do
    base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
    output_file=${base}_uniquely_mapped.bam
    echo "Filtering unique reads for $file"
    samtools view -@ 16 -q 255 -b $file > $output_file
done

samtools index -M *_uniquely_mapped.bam

echo "stringtie TPM estimates"
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
python3 03_get_TPM_table.py -t .
```
Outputs from 03_get_TPM_table.py can be found in /outputs/

Average_TPM.R simply averages the TPM across all the GRC genes in this data set 
