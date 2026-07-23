# Gene-rich fungus gnat germline-restricted chromosomes are largely transcriptionally quiescent

#### _Riccardo G. Kyriacou<sup>1,†</sup>, Marion Herbette<sup>2</sup>, Robert B. Baird<sup>3</sup>, Katy M. Monteith<sup>1</sup>, Yukiko M. Yamashita<sup>3,4</sup>, Laura Ross<sup>1,†,\*</sup>, Christina N. Hodson<sup>5,\*</sup>_

##### <sup>1</sup> Institute of Ecology and Evolution, University of Edinburgh, UK
##### <sup>2</sup> Laboratory of Biology and Modelling of the Cell, Ecole Normale Supérieure de Lyon, CNRS, Inserm, Université Claude Bernard Lyon 1, Lyon, France
##### <sup>3</sup> Whitehead Institute for Biomedical Science and Howard Hughes Medical Institute, Cambridge, Massachusetts, USA
##### <sup>4</sup> Department of Biology, Massachusetts Institute of Technology, Cambridge, MA, USA
##### <sup>5</sup> University College London, UCL Department of Genetics, Evolution & Environment, UK

###### \* Authors contributed equally
###### † Corresponding authors: R.G.Kyriacou@sms.ed.ac.uk, Laura.Ross@ed.ac.uk

---

> **Preprint:** [bioRxiv 2025.12.11.693641](https://www.biorxiv.org/content/10.64898/2025.12.11.693641v1)  
> **Note:** This repository accompanies the **revised manuscript**. Code associated with the original preprint is archived separately at [GRC_transcription](https://github.com/RiccardoKyriacou/GRC_transcription).

---

This repository provides scripts, analyses, and outputs for RNA-seq mapping, intergenic expression threshold estimation, GRC gene expression analysis, horizontal gene transfer (HGT) analysis, and re-analysis of published embryo RNA-seq data. Directories are ordered to follow the structure of the manuscript. README files within each subdirectory provide further usage information.

---

### Dependencies

#### Command-line tools

| Tool | Version | Used in |
|------|---------|---------|
| [fastp](https://github.com/OpenGene/fastp) | 0.24.0 | `01_RNAseq_Mapping` |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.12.1 | `01_RNAseq_Mapping` |
| [STAR](https://github.com/alexdobin/STAR) | 2.7.11b | `01_RNAseq_Mapping`, `05_Urban_et_al_re-analysis` |
| [StringTie](https://github.com/gpertea/stringtie) | 2.2.3 | `01_RNAseq_Mapping`, `02_Intergenic_mapping`, `05_Urban_et_al_re-analysis` |
| [samtools](https://www.htslib.org/) | 1.21 | `01_RNAseq_Mapping` |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | 2.16.0 | `03_Expressed_GRC_genes`, `04_GRC2_HGT` |
| [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/) | 5.76-107.0 | `03_Expressed_GRC_genes` |
| [minimap2](https://github.com/lh3/minimap2) | 2.28-r1209 | `04_GRC2_HGT` |
| [FASTGA](https://github.com/thegenemyers/FASTGA) | 1.1 | `04_GRC2_HGT` |

#### R packages

| Package | Used in |
|---------|---------|
| [mclust](https://cran.r-project.org/web/packages/mclust/index.html) | `02_Intergenic_mapping` |
| [ggplot2](https://ggplot2.tidyverse.org/) | `figures`, `supplementary_materials` |
| [patchwork](https://patchwork.data-imaginist.com/) | `figures`, `supplementary_materials` |

#### Python
All Python scripts are written for **Python 3**.

---

### Directory Summary

| Directory | Description |
|-----------|-------------|
| `01_RNAseq_Mapping` | Scripts to trim, quality-check, align RNA-seq reads to the reference _B. coprophila_ genome (GCA_965233685.1) using STAR, filter multi-mapped reads with samtools, and calculate TPM values with StringTie. Covers 42 libraries across four developmental stages (0–4h embryo, 4–8h embryo, late larva/early pupa, adult) in both germline and somatic libraries. |
| `02_Intergenic_mapping` | Scripts to generate intergenic GFF files (regions ≥500 bp from any annotated gene, ≥1 kb in length), map reads to intergenic regions, quantify TPM values, and apply a Gaussian mixture model (mclust) to estimate the background mismapping threshold (MITT ~0.55 TPM). |
| `03_Expressed_GRC_genes` | Scripts to identify confidently expressed GRC-linked genes above the active expression threshold (MATT), BLAST GRC transcripts against the core genome to flag high-similarity paralogues, and perform protein-level homology analysis against NCBI NR (BLASTp), Repbase, and InterProScan. Also includes Alien Index calculations to infer cecidomyiid vs. sciarid ancestry of expressed loci. |
| `04_GRC2_HGT` | Scripts to re-map PacBio long reads to the _B. coprophila_ assembly to validate the ~290 kb Rickettsiaceae-derived HGT region on GRC2, calculate GC content across the region, align GRC2 against the co-assembled Rickettsiaceae genome (FASTGA), and visualise synteny. |
| `05_Urban_et_al_re-analysis` | Scripts to download, align, and quantify pooled embryo RNA-seq data from Urban et al. (2021; BioProject PRJNA291918) spanning 2h–2 days post-fertilisation, to assess GRC transcription across and after zygotic genome activation. |
| `figures` | R scripts and outputs used to generate main paper Figures 2, 3, and 4. |
| `supplementary_materials` | R scripts and outputs used to generate supplementary figures and perform associated statistics. |

---

### Directory Structure

```text
Bcop_GRC_quiescence/
├── 01_RNAseq_Mapping/
│   ├── 01_STAR_TPM_Nmax2.sh
│   ├── 02_get_TPM_values.py
│   └── outputs/
│       ├── combined_TPM_only.tsv
│       ├── gene_overlap_0-4h_TPM.tsv
│       ├── gene_overlap_4-8h_TPM.tsv
│       ├── gene_overlap_adult_TPM.tsv
│       └── gene_overlap_late-larva-early-pupa_TPM.tsv
├── 02_Intergenic_mapping/
│   ├── 01_get_intergenic_GFF3.py
│   ├── 02_intergenic_TPM.sh
│   ├── 03_get_intergenic_TPM.py
│   ├── 04_intergenic_TPM_deconvolution.R
│   └── outputs/
│       ├── bcop_core_GRC.intergenic.gff3
│       ├── bcop_core_GRC.intergenic.gtf
│       └── combined_intergenic_TPM.tsv
├── 03_Expressed_GRC_genes/
│   ├── 01_get_expressed_genes.py
│   ├── 02_BLAST_GRCgenes.sh
│   ├── 03_get_BLAST_table.py
│   ├── 04_get_GRC_proteins.py
│   ├── 05_BLAST_GRC_proteins.sh
│   ├── 06_get_interpro_summary.py
│   ├── 07_repbase_tBLASTn_GRCproteins.sh
│   └── outputs/
│       ├── GRC_BLAST_table.tsv
│       ├── GRC_gene_expression.tsv
│       ├── GRC_gene_homology_BLAST/
│       │   ├── GRC_transcripts_BLAST_output.tsv
│       │   ├── GRCtranscript_repbase_output.tsv
│       │   └── iprscan5_GRC_genes.tsv
│       ├── GRC_v_Core_BLAST/
│       │   ├── GRC_v_Core_gene_BLAST_output.tsv
│       │   └── GRC_v_Core_genome_BLAST_output.tsv
│       └── fasta_files/
│           ├── GRC_genes.nucl.fasta
│           └── GRC_transcripts.fasta
├── 04_GRC2_HGT/
│   ├── 01_align_raw_reads_to_assembly.sh
│   ├── 02_BLAST_flanking_genes.sh
│   ├── 03_align_Rickettsia_assembly.sh
│   ├── 04_FastGA.sh
│   ├── 05_HGTregion_alignment.R
│   ├── 06_get_GC_content.py
│   ├── 07_sliding_window_GC.R
│   ├── get_kingdom.py
│   ├── data/
│   │   ├── 290kb_HGT_only.fasta
│   │   ├── 651kb_HGT_region.fasta
│   │   ├── 651kb_length.tsv
│   │   ├── GC_290kb_HGT_region_only_w2000_n1000.tsv
│   │   ├── GC_651kb_full_HGT_region_w2000_n1000.tsv
│   │   ├── GC_SUPER_GRC2_w200000_n100000.tsv
│   │   ├── Rickettsiaceae_contig_sizes.tsv
│   │   └── idBraCopr2.1.chrom_sizes.tsv
│   └── outputs/
│       └── grc2_vs_rickettsia_1to1.1aln.paf
├── 05_Urban_et_al_re-analysis/
│   ├── 01_download_pooled_embryo.sh
│   ├── 02_STAR_Stringtie.sh
│   ├── 03_get_TPM_table.py
│   ├── Average_TPM.R
│   └── outputs/
│       ├── combined_TPM_Urban_embryo.tsv
│       └── combined_TPM_Urban_embryo_GRC.tsv
├── figures/
│   ├── Figure_02.R
│   ├── Figure_03.R
│   ├── Figure_04.R
│   └── outputs/
│       ├── Figure_02.svg
│       ├── Figure_03.svg
│       └── Figure_04.svg
├── supplementary_materials/
│   ├── Figure_S2.R
│   ├── Figure_S3.R
│   └── outputs/
│       ├── Figure_S2.svg
│       └── Figure_S3.svg
└── README.md
