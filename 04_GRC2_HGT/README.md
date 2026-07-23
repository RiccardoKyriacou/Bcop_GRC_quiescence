# 4. HGT analysis

RNA-seq analysis revealed two genes, g19161 and g19121, with clear bacterial origin.These scripts analyse the raw read libraries and identify a ~290kb bacterial region in GRC2, which may have evolved through HGT with endosymbiotic _Rickettsia_ 

First, using 01_align_raw_reads_to_assembly.sh, we run minimap2 to index and then re-map the raw, HiFi long-reads back to the _B coprophila_ genome.
```
minimap2 -d idBraCopr2.1.primary.masked.mmi idBraCopr2.1.primary.masked.fa
minimap2 -ax map-pb idBraCopr2.1.primary.masked.mmi ERR12736861.fastq.gz > ERR12736861.sam

# Convert to sorted BAM
samtools view -@ 16 -bS ERR12736861.sam | samtools sort -@ 16 -o ERR12736861.sorted.bam
samtools index ERR12736861.sorted.bam
rm ERR12736861.sam
```
We next use 02_BLAST_flanking_genes.sh to BLAST a large subset of gene (g19064-g19168), which flank this HGT region, to better understand which genes are eukaryotic and which are bacterial 
```
blastp \
-query /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs/flanking_gene_transcripts.fasta \
-db /mnt/loki/db/ncbi_nr/nr \
-outfmt "6 std sscinames staxids stitle" \
-num_threads 16 \
-max_target_seqs 5 \
-out /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs/flanking_genes_BLAST_output.tsv
```
We also use 03_align_Rickettsia_assembly.sh to align the _Rickettsia_ genome to the full _B. coprophila_ reference genome using minimiap2
```
minimap2 -ax asm5 idBraCopr2.1.primary.masked.fa Rickettsiaceae.finalassembly.fa > Bcop_Rickettsia.sam

echo "Converting to BAM, sorting, and indexing..."
samtools view -bS Bcop_Rickettsia.sam | samtools sort -o Bcop_Rickettsia.sorted.bam
samtools index Bcop_Rickettsia.sorted.bam
```
04_FastGA.sh is also used to run FASTGA to generate a PAF alignment file between the two genomes (see /outputs/grc2_vs_rickettsia_1to1.1aln.paf)
```
# Run FastGA (produce ALN file for plotting)
FastGA \
  Rickettsiaceae.finalassembly.fa \
  idBraCopr2.1.primary.masked.fa \
  -1:grc2_vs_rickettsia.1aln \
  -T16 -vk

# Force 1 to 1 alignment
ALNchain -v -ogrc2_vs_rickettsia_1to1.1aln grc2_vs_rickettsia.1aln
 
# Convert to PAF
ALNtoPAF grc2_vs_rickettsia_1to1.1aln > grc2_vs_rickettsia_1to1.1aln.paf
```

We also look at GC content across the HGT region and GRC2 as a whole. To do this we use the same python script, but increase the parameters of the sliding window. The script requires a genome (or any DNA sequence) (-g), a window size in bp (-w) and a step size in bp (-n).

To get the GC content across the HGT region we run the script first across a 651kb region containing the HGT region (see /data/651kb_HGT_region.fasta) and then only across the HGT region itself (/data/290kb_HGT_only.fasta)
```
05_get_GC_content.py -g /data/651kb_HGT_region.fasta - w 2000 -n 1000
05_get_GC_content.py -g /data/290kb_HGT_only.fasta - w 2000 -n 1000
```
We also run the script across the whole of GRC2, but increase the window and step size 
```
05_get_GC_content.py -g idBraCopr2.1.primary.masked.fa -w 200000 -n 100000
```
GC outputs are all availible in /outputs/GC_content
