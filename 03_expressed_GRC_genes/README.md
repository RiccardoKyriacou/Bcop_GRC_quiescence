# 3. Defining confidently expressed GRC-linked genes and downstream analysis

## Overview
This directory contains scripts to identify confidently expressed GRC-linked genes in 
_B. coprophila_, characterise their similarity to core-chromosome paralogues, and perform 
downstream homology analysis to investigate potential function and ancestry. This includes 
nucleotide and protein-level BLAST searches, InterProScan domain prediction, Repbase TE 
homology searches, and Alien Index calculations to infer cecidomyiid vs. sciarid gene ancestry.

---

## Expression filtering criteria

For a GRC-linked gene to be considered **confidently expressed** in this study, it must 
satisfy **all** of the following conditions:

- **TPM > ~0.55** (threshold derived from `02_intergenic_TPM`) in at least **2/3 of 
replicate libraries** per sex per developmental stage
- **No significant expression in matched somatic libraries** — as GRCs are absent from 
somatic cells, any apparent somatic expression reflects core chromosome reads mismapping 
to paralogous GRC loci. These same mismapped reads could falsely inflate expression 
estimates in germline libraries, so affected genes are excluded
- **No highly similar core-chromosome paralogue** — defined as a BLASTn hit with 
% identity × % coverage > 70, as reads from highly similar core paralogues could 
misalign to GRC loci and generate false-positive expression signals

---

## Scripts

### Step 1 — Classify expressed GRC-linked genes

Loading `TPM_genes.tsv` (`Bcop_GRC_quiescence/02_intergenic_TPM/outputs/TPM_genes.tsv`) 
into R, we run:

```
01_classify_GRC_expressed_genes.R
```

This script applies the expression filtering criteria above across all developmental 
stages and both sexes. Genes must exceed the TPM threshold in at least 2/3 of replicates 
per sex per stage. It generates the following output tables saved to 
`outputs/GRC_classification_outputs/`:

| Output file | Description |
|-------------|-------------|
| `Expressed_GRC_genes.tsv` | Full table of all GRC-linked genes meeting the TPM threshold in at least one developmental stage/sex |
| `Expressed_GRC_genes_summary.tsv` | Summarised table of expressed GRC genes, collapsed across stages |
| `Corrected_expressed_GRC_genes.tsv` | Final filtered table with genes showing significant somatic library mismapping removed |

---

### Step 2 — Extract GRC and core gene nucleotide sequences

Using `Corrected_expressed_GRC_genes.tsv` to define our gene set, we extract nucleotide 
sequences for both the confidently expressed GRC-linked genes and all annotated core 
chromosome genes in _B. coprophila_. Core gene sequences are needed as the BLAST database 
in Step 3.

```bash
python 02_get_GRC_core_fastas.py \
    -g /path/to/Annotations/idBraCopr2.1.primary.masked.fa \
    -a /path/to/Annotations/bcop_core_GRC.gff3 \
    -s outputs/GRC_classification_outputs/Corrected_expressed_GRC_genes.tsv
```

| Argument | Description |
|----------|-------------|
| `-g` | Path to masked reference genome FASTA |
| `-a` | Path to GFF3 annotation |
| `-s` | Path to `Corrected_expressed_GRC_genes.tsv` from Step 1 |

| Output file | Description |
|-------------|-------------|
| `outputs/fasta_files/GRC_genes_to_BLAST.fasta` | Full gene body nucleotide sequences (including introns) for confidently expressed GRC-linked genes |
| `outputs/fasta_files/core_genes.fasta` | Full gene body nucleotide sequences for all annotated core chromosome genes |

> **Note:** This script extracts full gene body sequences (including introns) rather than 
> spliced CDS. These are used for nucleotide-level BLAST comparisons in Steps 3–4 only. 
> Spliced CDS sequences are generated separately in Step 5 for protein-level analysis.

---

### Step 3 — BLASTn GRC genes against core chromosomes

We perform BLASTn searches of expressed GRC gene sequences against both the annotated 
core chromosome genes and the full masked core genome. This identifies GRC-linked genes 
with highly similar core chromosome paralogues, which would reduce confidence in their 
expression due to potential mismapping of core reads to GRC loci.

```bash
03_BLAST_GRCgenes.sh
```

**Step 3a** — BLASTn against annotated core gene sequences:

```bash
# Make BLAST database of core chromosome genes
makeblastdb \
    -in core_genes.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out core_genes_DB

# BLASTn GRC genes against core genes
blastn \
    -query GRC_genes_to_BLAST.fasta \
    -db core_genes_DB \
    -out GRC_v_Core_gene_BLAST_output.tsv \
    -outfmt '6 std qlen slen qseq sseq'
```

**Step 3b** — BLASTn against the masked core genome (to capture unannotated paralogues):

A custom Python script (`mask_genome.py`) first masks all GRC scaffolds in the reference 
genome, ensuring only core chromosome sequence is searched against:

```bash
# Mask GRC scaffolds from the genome
python3 mask_genome.py \
    -g idBraCopr2.1.primary.masked.fa \
    -o masked_GRC_genome.fasta \
    -m GRC

# Make BLAST database of masked core genome
makeblastdb \
    -in masked_GRC_genome.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out masked_GRC_genome_DB

# BLASTn GRC genes against masked core genome
blastn \
    -query GRC_genes_to_BLAST.fasta \
    -db masked_GRC_genome_DB \
    -out GRC_v_Core_genome_BLAST_output.tsv \
    -outfmt '6 std qlen slen qseq sseq'
```

Outputs are saved to `outputs/GRC_v_Core_BLAST/`.

---

### Step 4 — Generate BLAST summary table and flag high-similarity core paralogues

BLAST results from Step 3 are passed through a custom Python script to generate a summary 
table of expressed GRC-linked genes alongside their best-hit core chromosome paralogue. 
A similarity score (% identity × % coverage) is calculated for each hit. Genes with a 
score > 70 are flagged as having a highly similar core paralogue and are excluded from 
the final confident set, as core chromosome reads could plausibly misalign to these 
GRC loci.

```bash
python 04_get_BLAST_table.py \
    -t outputs/GRC_classification_outputs/Corrected_expressed_GRC_genes.tsv \
    -c outputs/GRC_v_Core_BLAST/GRC_v_Core_gene_BLAST_output.tsv \
    -g outputs/GRC_v_Core_BLAST/GRC_v_Core_genome_BLAST_output.tsv
```

| Argument | Description |
|----------|-------------|
| `-t` | Path to `Corrected_expressed_GRC_genes.tsv` from Step 1 |
| `-c` | Path to GRC vs. core gene BLAST output from Step 3a |
| `-g` | Path to GRC vs. core genome BLAST output from Step 3b |

| Output | Description |
|--------|-------------|
| `outputs/GRC_BLAST_table.tsv` | Table of expressed GRC genes with best-hit core paralogue similarity scores |

---

### Step 5 — Extract and translate spliced CDS sequences to protein

Analysis now focuses on the confidently expressed GRC-linked genes. To perform 
protein-level homology analysis, we extract properly spliced CDS sequences and translate 
them to amino acid sequences.

Unlike Step 2 which extracted unspliced gene body sequences, this script:
- Parses only `CDS` features from the GFF annotation (excluding introns)
- Groups CDS segments by transcript ID
- Splices segments together in genomic coordinate order
- Applies **phase correction** by trimming the spliced sequence by the phase offset of 
the biologically first CDS, ensuring the sequence is in the correct reading frame
- Handles strand orientation correctly for both `+` and `-` strand transcripts
- Translates the in-frame spliced CDS to amino acid sequences

```bash
python 05_get_GRC_proteins.py \
    -t outputs/GRC_classification_outputs/Corrected_expressed_GRC_genes.tsv \
    -g /path/to/Annotations/idBraCopr2.1.primary.masked.fa \
    -a /path/to/Annotations/bcop_core_GRC.gff3
```

| Argument | Description |
|----------|-------------|
| `-t` | Path to `Corrected_expressed_GRC_genes.tsv` from Step 1 |
| `-g` | Path to masked reference genome FASTA |
| `-a` | Path to GFF3 annotation |

| Output | Description |
|--------|-------------|
| `outputs/fasta_files/GRC_transcripts.fasta` | Translated amino acid sequences for expressed GRC-linked transcripts |
| `outputs/fasta_files/GRC_genes.nucl.fasta` | Spliced in-frame nucleotide CDS sequences for expressed GRC-linked transcripts |

---

### Step 6 — BLASTp against NCBI non-redundant protein database

Protein sequences are queried against the NCBI non-redundant (NR) protein database 
(BLASTp v2.16.0) to identify similarities to known proteins across all taxa:

```bash
06_BLAST_GRC_proteins.sh
```

| Output | Description |
|--------|-------------|
| `outputs/GRC_gene_homology_BLAST/GRC_transcripts_BLAST_output.tsv` | BLASTp hits for expressed GRC protein sequences against NCBI NR |

---

### Step 7 — tBLASTn against Repbase to assess TE homology

To assess whether expressed GRC-linked loci show homology to transposable elements, 
protein sequences are queried against the Repbase database of repetitive elements 
(Bao et al. 2015) using tBLASTn (protein query vs. nucleotide database):

```bash
07_repbase_tBLASTn_GRCproteins.sh
```

| Output | Description |
|--------|-------------|
| `outputs/GRC_gene_homology_BLAST/GRCtranscript_repbase_output.tsv` | tBLASTn hits for expressed GRC protein sequences against Repbase |

---

### Step 8 — InterProScan domain prediction

All expressed GRC protein sequences are submitted to InterProScan (v5.76-107.0) to 
identify conserved protein domains and assign protein family classifications. Results 
were obtained via the [EBI InterProScan web service](https://www.ebi.ac.uk/interpro/search/sequence/) 
and parsed using a custom Python script:

```bash
python 08_get_interpro_summary.py
```

| Output | Description |
|--------|-------------|
| `outputs/GRC_gene_homology_BLAST/iprscan5_GRC_genes.tsv` | InterProScan domain predictions and protein family classifications for expressed GRC protein sequences |

---

## Outputs summary

```
outputs/
├── GRC_classification_outputs/
│   ├── Expressed_GRC_genes.tsv
│   ├── Expressed_GRC_genes_summary.tsv
│   └── Corrected_expressed_GRC_genes.tsv
├── GRC_BLAST_table.tsv
├── GRC_gene_expression.tsv
├── fasta_files/
│   ├── GRC_genes_to_BLAST.fasta
│   ├── core_genes.fasta
│   ├── GRC_genes.nucl.fasta
│   └── GRC_transcripts.fasta
├── GRC_v_Core_BLAST/
│   ├── GRC_v_Core_gene_BLAST_output.tsv
│   └── GRC_v_Core_genome_BLAST_output.tsv
└── GRC_gene_homology_BLAST/
    ├── GRC_transcripts_BLAST_output.tsv
    ├── GRCtranscript_repbase_output.tsv
    └── iprscan5_GRC_genes.tsv
```
```

---

### Things to double-check before committing

| Item | Note |
|------|------|
| Script numbering (Steps 5–8) | Confirm the actual filenames in the repo match — these were inferred from the original directory listing and paper methods |
| `GRC_gene_expression.tsv` | This output appears in the directory but is not clearly generated by any of the documented steps — confirm which script produces it |
| Cross-species BLASTn | The paper describes BLASTn searches against _L. ingenua_ and _B. impatiens_ GRCs (Figure 3C) — confirm whether this script exists in the directory and add a step if so |
| Three-way BLASTp / Alien Index | The paper describes a three-way BLASTp comparison against _A. aphidimyza_, _B. coprophila_ core, and _A. gambiae_ proteomes to calculate the Alien Index — confirm whether this is in a separate script and add a step if so |
