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
> **Note:** This repository accompanies a  **manuscript currently under revision**. For code associated with the original preprint, see [GRC_transcription](https://github.com/RiccardoKyriacou/GRC_transcription).

---

Directory provides scripts, analyses, and outputs for RNA-seq mapping, intergenic mapping, GRC gene expression, and horizontal gene transfer analysis, ordered accordingly to the manuscript.  
README files provided for each subdirectory providing usage information for scripts.

---

### Dependencies

| Tool | Version | Used in |
|------|---------|---------|
| STAR | x.x.x | `01_RNAseq_Mapping`, `05_Urban_et_al_re-analysis` |
| StringTie | x.x.x | `01_RNAseq_Mapping`, `05_Urban_et_al_re-analysis` |
| BLAST+ | x.x.x | `03_Expressed_GRC_genes`, `04_GRC2_HGT` |
| Python | 3.x | Multiple |
| R | 4.x.x | Multiple |

---

### Directory Summary

| Directory | Description |
|-----------|-------------|
| `01_RNAseq_Mapping` | Scripts to align RNA-seq reads to the reference _B. coprophila_ genome using STAR and calculate TPM values with StringTie. |
| `02_Intergenic_mapping` | Scripts to generate intergenic GFF files, map reads to intergenic regions, and calculate intergenic TPMs. |
| `03_Expressed_GRC_genes` | Scripts to identify expressed GRC-linked genes, BLAST GRC genes against core genome, and perform homology analysis against NCBI NR proteins and RepBase. |
| `04_HGT` | Scripts to investigate putative horizontal gene transfer (HGT) regions in GRC2. |
| `05_Urban_et_al_re-analysis` | Scripts to re-analyse pooled embryo RNA-seq data from Urban et al. (2021). |
| `figures` | Raw figures and code used to generate main paper Figures 2, 3 and 4. |
| `supplementary_materials` | Raw figures and code used to generate supplementary figures and perform statistics. |

---

### Data Availability
Raw sequencing reads and genome assembly are deposited at [NCBI/ENA — accession numbers to be added].

---

### Citation
Kyriacou et al. (2026) Gene-rich fungus gnat germline-restricted chromosomes are largely transcriptionally quiescent. *[Journal]*. DOI: [to be added upon publication]
