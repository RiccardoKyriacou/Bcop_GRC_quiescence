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
> **Note:** This repository accompanies **a manuscript currently under revision**. Code associated with the original preprint is archived separately at [GRC_transcription](https://github.com/RiccardoKyriacou/GRC_transcription).

---

This repository provides scripts, analyses, and outputs for RNA-seq mapping, intergenic expression threshold estimation, GRC gene expression analysis, horizontal gene transfer (HGT) analysis, and re-analysis of published embryo RNA-seq data. Directories are ordered to follow the structure of the manuscript. README files within each subdirectory provide further usage information.

---

### Dependencies

#### Command-line Tools

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

#### R Packages

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
| `02_Intergenic_mapping` | Scripts to generate intergenic GFF files (regions ≥500 bp from any annotated gene, ≥1 kb in length), quantify intergenic TPM values, and apply a Gaussian mixture model (mclust) to estimate the background mismapping threshold (~0.55 TPM). |
| `03_Expressed_GRC_genes` | Scripts to identify confidently expressed GRC-linked genes above the active expression threshold (MATT), BLAST GRC transcripts against the core genome to flag high-similarity paralogues, and perform protein-level homology analysis against NCBI NR (BLASTp), Repbase, and InterProScan. Also includes Alien Index calculations to infer cecidomyiid vs. sciarid ancestry of expressed loci. |
| `04_GRC2_HGT` | Scripts to re-map PacBio long reads to the _B. coprophila_ assembly to validate the ~290 kb Rickettsiaceae-derived HGT region on GRC2, calculate GC content across the region in a sliding window, align GRC2 against the co-assembled Rickettsiaceae genome using FASTGA, and visualise synteny. |
| `05_Urban_et_al_re-analysis` | Scripts to download, align, and quantify pooled embryo RNA-seq data from Urban et al. (2021; BioProject PRJNA291918) spanning 2h–2 days post-fertilisation, to assess GRC transcription across and after zygotic genome activation. |
| `figures` | R scripts and outputs used to generate main paper Figures 2, 3, and 4. |
| `supplementary_materials` | R scripts and outputs used to generate supplementary figures and perform associated statistics. |

---

### Data Availability

| Data | Accession |
|------|-----------|
| Reference _B. coprophila_ genome assembly | [GCA_965233685.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_965233685.1/); SRA: ERS15411730 |
| Adult RNA-seq libraries (Baird et al. 2025) | BioProject [PRJNA1109384](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1109384) |
| Embryo RNA-seq libraries (this study) | BioProject [PRJNA1220056](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1220056) |
| Late larval/early pupal RNA-seq (this study) | Available upon publication |
| _Rickettsiaceae_ endosymbiont genome | Available upon publication |
| Urban et al. (2021) pooled embryo RNA-seq | BioProject [PRJNA291918](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA291918) |

---

### Citation

> Kyriacou RG, Herbette M, Baird RB, Monteith KM, Yamashita YM, Ross L, Hodson CN. (2026). Gene-rich fungus gnat germline-restricted chromosomes are largely transcriptionally quiescent. *[Journal]*. DOI: [to be added upon publication]

---

### License

This repository is licensed under the [MIT License](LICENSE).
