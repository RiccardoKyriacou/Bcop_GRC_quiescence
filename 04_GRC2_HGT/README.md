# 4. HGT analysis

## Overview

RNA-seq analysis of GRC-linked genes revealed two loci, g19161 and g19121, with clear 
bacterial origin — g19161 showing homology to the chaperone cofactor GrpE, and g19121 
to various bacterial hypothetical proteins, with closest hits in Rickettsiaceae. These 
scripts investigate whether these genes represent a genuine horizontal gene transfer 
(HGT) event from the endosymbiotic Rickettsiaceae known to inhabit _B. coprophila_ 
germline cells, or erroneous assembly of a contaminating bacterial contig into the 
reference genome.

---

## Scripts

### Step 1 — Re-map HiFi long reads to the _B. coprophila_ genome

Using minimap2 (v2.28-r1209), we index the reference genome and re-map the original 
PacBio HiFi long reads used to assemble it

```bash
01_align_raw_reads_to_assembly.sh
```

```bash
# Index the reference genome
minimap2 -d idBraCopr2.1.primary.masked.mmi idBraCopr2.1.primary.masked.fa

# Map HiFi long reads
minimap2 -ax map-pb idBraCopr2.1.primary.masked.mmi ERR12736861.fastq.gz > ERR12736861.sam

# Convert to sorted BAM
samtools view -@ 16 -bS ERR12736861.sam | samtools sort -@ 16 -o ERR12736861.sorted.bam
samtools index ERR12736861.sorted.bam
rm ERR12736861.sam
```

| Output | Description |
|--------|-------------|
| `ERR12736861.sorted.bam` | Sorted BAM file of HiFi reads aligned to the _B. coprophila_ reference genome |
| `ERR12736861.sorted.bam.bai` | BAM index for visualisation in IGV |

---

### Step 2 — BLASTp flanking genes to classify as eukaryotic or bacterial

To characterise the ~290 kb region and its immediate flanking sequence, we BLASTp 
predicted protein sequences for a broader set of genes spanning the region (g19064–g19168) 
against the NCBI NR protein database. This allows each gene to be classified as 
eukaryotic or bacterial in origin, and defines the boundaries of the HGT insertion.

```bash
02_BLAST_flanking_genes.sh
```

```bash
blastp \
    -query flanking_gene_transcripts.fasta \
    -db /path/to/ncbi_nr/nr \
    -outfmt "6 std sscinames staxids stitle" \
    -num_threads 16 \
    -max_target_seqs 5 \
    -out flanking_genes_BLAST_output.tsv
```

| Output | Description |
|--------|-------------|
| `outputs/flanking_genes_BLAST_output.tsv` | BLASTp hits for genes g19064–g19168 against NCBI NR, including taxonomic annotations |

---

### Step 3 — Align Rickettsiaceae genome to _B. coprophila_ reference

Using minimap2 with the `asm5` preset (suitable for divergent genome-to-genome 
alignments), we align the co-assembled Rickettsiaceae endosymbiont genome against 
the full _B. coprophila_ reference to identify regions of synteny between the 
endosymbiont and the putative GRC2 HGT region:

```bash
03_align_Rickettsia_assembly.sh
```

```bash
minimap2 -ax asm5 idBraCopr2.1.primary.masked.fa Rickettsiaceae.finalassembly.fa > Bcop_Rickettsia.sam

samtools view -bS Bcop_Rickettsia.sam | samtools sort -o Bcop_Rickettsia.sorted.bam
samtools index Bcop_Rickettsia.sorted.bam
```

| Output | Description |
|--------|-------------|
| `Bcop_Rickettsia.sorted.bam` | Sorted BAM of Rickettsiaceae genome aligned to _B. coprophila_ reference |
| `Bcop_Rickettsia.sorted.bam.bai` | BAM index file |

---

### Step 4 — Generate genome-to-genome alignment with FASTGA

FASTGA (v1.1) is used to generate a high-quality pairwise alignment between the 
Rickettsiaceae genome and the _B. coprophila_ reference genome. A 1-to-1 chain 
alignment is forced to remove redundant overlapping hits, before converting to PAF 
format for downstream synteny visualisation and divergence estimation in Steps 5–6.

```bash
04_FastGA.sh
```

```bash
# Run FASTGA
FastGA \
    Rickettsiaceae.finalassembly.fa \
    idBraCopr2.1.primary.masked.fa \
    -1:grc2_vs_rickettsia.1aln \
    -T16 -vk

# Force 1-to-1 alignment
ALNchain -v -ogrc2_vs_rickettsia_1to1.1aln grc2_vs_rickettsia.1aln

# Convert to PAF
ALNtoPAF grc2_vs_rickettsia_1to1.1aln > grc2_vs_rickettsia_1to1.1aln.paf
```

| Output | Description |
|--------|-------------|
| `outputs/grc2_vs_rickettsia_1to1.1aln.paf` | 1-to-1 PAF alignment between Rickettsiaceae genome and _B. coprophila_ GRC2, used for synteny visualisation and divergence estimation |

---

### Step 5 — Calculate GC content across the HGT region and GRC2

A custom Python script calculates GC content in a sliding window across a given DNA 
sequence. The bacterial-derived HGT region is expected to show a visible dip in GC 
content relative to the surrounding GRC2 sequence, consistent with the typically lower 
GC content of bacterial genomes. The script is run at three different scales to 
capture this pattern.

| Argument | Description |
|----------|-------------|
| `-g` | Path to input FASTA file |
| `-w` | Sliding window size (bp) |
| `-n` | Step size (bp) |

**Across the 651 kb region containing the HGT insertion:**

```bash
python 05_get_GC_content.py -g data/651kb_HGT_region.fasta -w 2000 -n 1000
```

**Across the 290 kb HGT region only:**

```bash
python 05_get_GC_content.py -g data/290kb_HGT_only.fasta -w 2000 -n 1000
```

**Across the whole of GRC2** (larger window to capture the broader GC landscape):

```bash
python 05_get_GC_content.py -g idBraCopr2.1.primary.masked.fa -w 200000 -n 100000
```

| Output | Description |
|--------|-------------|
| `data/GC_651kb_full_HGT_region_w2000_n1000.tsv` | GC content across the 651 kb region in 2 kb windows |
| `data/GC_290kb_HGT_region_only_w2000_n1000.tsv` | GC content across the 290 kb HGT region only in 2 kb windows |
| `data/GC_SUPER_GRC2_w200000_n100000.tsv` | GC content across GRC2 in 200 kb windows |

---

### Step 6 — Estimate age of the HGT event from alignment divergence

To estimate when the Rickettsiaceae-derived region was integrated into GRC2, we 
calculate nucleotide divergence between the _B. coprophila_ HGT sequence and its 
corresponding region in the Rickettsiaceae genome using the 1-to-1 PAF alignment 
generated in Step 4.

Divergence (D) is calculated as the total counted nucleotide differences across all 
alignment blocks (extracted from the `df:i` tag in the FASTGA PAF output) divided by 
the total alignment length:

$$D = \\frac{\\sum \\text{differences}}{\\sum \\text{alignment length}}$$

An estimated divergence time (T) in generations is then derived assuming neutral 
evolution across the entire region:

$$T = \\frac{D}{\\mu_{\\text{fly}} + \\mu_{\\text{bacteria}}}$$

Since no mutation rate estimates exist for _B. coprophila_ or its Rickettsiaceae 
endosymbiont, we substitute reported neutral mutation rates for _Drosophila melanogaster_ 
( $ \mu = 2.8 \times 10^{-9} $ per site per generation; Keightley et al. 2014) and the 
endosymbiotic bacterium _Teredinibacter turnerae_ 
( $ \mu = 1.1 \times 10^{-9} $ per site per generation; Senra et al. 2018).

```bash
python 06_get_divergence.py
```

> **Note:** The PAF file path is hardcoded in the script as 
> `grc2_vs_rickettsia_1to1.1aln.paf`. Ensure this file is present in the working 
> directory before running, or update the `PAF_FILE` variable at the top of the script.

The script prints divergence and estimated age in generations to the terminal. To 
convert to years, multiply by the reported _B. coprophila_ generation time of 
**24–40 days** (Baird et al. 2023):

$$T_{\text{years}} = T_{\text{generations}} \times \frac{\text{generation time (days)}}{365}$$



