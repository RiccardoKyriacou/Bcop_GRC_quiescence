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

### Step 1 — Extract GRC and core gene nucleotide sequences

Using `TPM_genes.tsv` to define our gene set, we extract nucleotide 
sequences for both the GRC-linked genes and all annotated core 
chromosome genes in _B. coprophila_.

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
| `outputs/fasta_files/GRC_genes_to_BLAST.fasta` | Full gene body nucleotide sequences (including introns) GRC-linked genes non-zero TPM |
| `outputs/fasta_files/core_genes.fasta` | Full gene body nucleotide sequences for all annotated core chromosome genes |

> **Note:** This script extracts full gene body sequences (including introns) rather than 
> spliced CDS. These are used for nucleotide-level BLAST comparisons in Steps 3–4 only. 
> Spliced CDS sequences are generated separately in Step 5 for protein-level analysis.

---

### Step 2 — BLASTn GRC genes against core chromosomes

We perform BLASTn searches of expressed GRC gene sequences against both the annotated 
core chromosome genes and the full masked core genome. This identifies GRC-linked genes 
with highly similar core chromosome paralogues, which would reduce confidence in their 
expression due to potential mismapping of core reads to GRC loci.

```bash
02_BLAST_GRCgenes.sh
```

**Step 2a** — BLASTn against annotated core gene sequences:

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

**Step 2b** — BLASTn against the masked core genome (to capture unannotated paralogues):

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

> ***Note*** to reduce file size, only results for confidently expressed GRC-linked genes shown (defined in the next step)

---

### Step 3 — Generate summary BLASTn table 

We can then generate a summary BLASTn file for GRC-linked genes with TPM > 0 using the following custom script:

```bash
python 03_get_BLAST_table.py \
    -t /path/to/TPM_genes.tsv \
    -c /path/to/core_gene_blast.tsv \
    -g /path/to/core_genome_blast.tsv \
    -o GRC_BLAST_table.tsv
```

| Argument | Description |
|----------|-------------|
| `-t` | Path to `Bcop_GRC_quiescence/02_intergenic_TPM/outputs/TPM_genes.tsv` |
| `-c` | Path to GRC genes v core gene BLASTn `/outputs/GRC_gene_BLAST_output` |
| `-g` | Path to GRC genes v core genome BLASTn `/outputs/GRC_genome_BLAST_output`|
| `-o` | Output file name `/outputs/GRC_BLAST_table.tsv`|

---

### Step 4 — Classify expressed GRC-linked genes

Loading `TPM_genes.tsv` (`Bcop_GRC_quiescence/02_intergenic_TPM/outputs/TPM_genes.tsv`) 
into R, we run:

```
04_classify_GRC_expressed_genes.R
```

This script also requires you top load in BLASTn summary table generated by `03_get_BLAST_table.py`

This script applies the expression filtering criteria above across all developmental 
stages and both sexes. Genes must exceed the TPM threshold in at least 2/3rds of replicates 
per sex per stage. They also must not have a core chromosome paralogue with a similarity that 
exceeds %identity * %coverage < 70. The script generates the following output tables saved to:
`outputs/expressed_GRC_gene_tables`:

| Output file | Description |
|-------------|-------------|
| `Expressed_GRC_genes.tsv` | Full table of all GRC-linked genes meeting the TPM threshold in at least one developmental stage/sex |
| `Expressed_GRC_genes_summary.tsv` | Summarised table of expressed GRC genes, collapsed across stages |
| `Corrected_expressed_GRC_genes.tsv` | Final filtered table with genes showing significant somatic library mismapping removed |

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

### Step 6 — BLASTp against NCBI non-redundant protein database and Repbase 

We now use a bash script to perform BLAST searches against both the NCBI (NR) protein database and Repbase

```bash
06_NCBI_repbase_BLAST.sh
```

**Step 6a** — BLASTp against NCBI non-redundant (NR) protein database:

Protein sequences are queried against the NCBI non-redundant (NR) protein database 
(BLASTp v2.16.0) to identify similarities to known proteins across all taxa:

```bash
echo "BLASTp for GRC proteins against NCBI nr"
blastp \
-query $WORKDIR/GRC_transcripts.fasta \
-db /mnt/loki/db/ncbi_nr/nr \
-outfmt "6 std sscinames staxids stitle" \
-evalue 0.05 \
-num_threads 16 \
-out $WORKDIR/GRC_transcripts_BLAST_output.tsv
```
**Step 6b** — tblastn against Repbase:

To assess whether expressed GRC-linked loci show homology to transposable elements, 
protein sequences are queried against the Repbase database of repetitive elements 
(Bao et al. 2015) using tBLASTn (protein query vs. nucleotide database):

```bash
echo "BLASTp for GRC proteins against Repbase"
tblastn \
-query $WORKDIR/GRC_transcripts.fasta \
-db /mnt/loki/db/Repbase25.09/Repbase25.09.fasta \
-outfmt "6 std" \
-evalue 1e-5 \
-num_threads 16 \
-out $WORKDIR/GRC_transcript_repbase_output.tsv
```

| Output | Description |
|--------|-------------|
| `outputs/GRC_transcript_blast/GRC_transcripts_BLAST_output.tsv` | BLASTp hits for expressed GRC protein sequences against NCBI NR |
| `outputs/GRC_transcript_blast/GRCtranscript_repbase_output.tsv` | tBLASTn hits for expressed GRC protein sequences against Repbase |

---

### Step 7 — InterProScan domain prediction

All expressed GRC protein sequences are submitted to InterProScan (v5.76-107.0) to 
identify conserved protein domains and assign protein family classifications. Results 
were obtained via the [EBI InterProScan web service](https://www.ebi.ac.uk/interpro/search/sequence/) 
and parsed using a custom Python script:

```bash
python 07_get_interpro_summary.py
```

| Output | Description |
|--------|-------------|
| `outputs/GRC_gene_homology_BLAST/iprscan5_GRC_genes.tsv` | InterProScan domain predictions and protein family classifications for expressed GRC protein sequences |

> ***Note*** Screenshots of Interpro domains included, these were used to produce Figure 3 Panel B `Bcop_GRC_quiescence/main_figures/Fig3`

---

### Step 8 — BLASTn of expressed GRC-linked genes against *B. impatiens* and *L. ingenua*

To assess cross-species conservation of confidently expressed GRC-linked loci, we perform 
BLASTn searches of expressed GRC nucleotide sequences against the full genome assemblies 
(GRCs and core chromosomes) of two other fungus gnat species with assembled GRC genomes: 
*Bradysia impatiens* and *Lycoriella ingenua* (Hodson et al. 2026).

```bash
08_Bimp_Ling_BLAST.sh
```

Nucleotide sequences from `outputs/fasta_files/GRC_genes.nucl.fasta` are queried against 
both species' assemblies. Only hits with >80% sequence identity and high query coverage 
are retained.

| Output | Description |
|--------|-------------|
| `outputs/Bimp_Ling_BLAST/` | BLASTn results for expressed GRC-linked loci against *B. impatiens* and *L. ingenua* assemblies |

> ***Note:*** Cross-species BLASTn results are used to produce Figure 3C. Results reveal 
> differing patterns of conservation: TE-like loci (e.g., Kolobok-like g13363, BEL-like 
> g17107) showed strong GRC homology in *B. impatiens*, the Transib-like locus (g596) was 
> detected on *B. impatiens* chromosome II, and the putative transposase-like locus (g7958) 
> showed homology to both the X chromosome and GRC1 of *L. ingenua*. In contrast, no GRC 
> homology was detected for insect-like genes in either species, though g18444 produced 
> strong hits to the core X chromosome of both *B. impatiens* and *L. ingenua*.

---

### Step 9 — Alien Index calculation to infer GRC gene ancestry

To investigate whether confidently expressed GRC-linked transcripts are more likely derived 
from an ancient hybridisation with a Cecidomyiidae ancestor or from subsequent duplication 
from the *B. coprophila* core genome, we calculate an Alien Index (AI) for each transcript. 
This step uses the BLASTp outputs generated in **Step 10** and should therefore be run 
after Step 10.

The Alien Index (Rancurel et al. 2017) is defined as:

$$AI = \ln(E_{\text{sciarid}} + 1 \times 10^{-200}) - \ln(E_{\text{cecidomyiid}} + 1 \times 10^{-200})$$

Higher AI values indicate greater similarity to cecidomyiid proteins (cecidomyiid-like 
ancestry), while lower values indicate greater similarity to the sciarid core proteome 
(likely derived from the *B. coprophila* core genome). Transcripts are classified as 
cecidomyiid-like (AI > 10), sciarid-like (AI < −10), or ambiguous.

> ***Note:*** This step requires the BLASTp outputs from **Step 10** 
> (`GRC_vs_SciaridCore.tsv`, `GRC_vs_Cecidomyiid.tsv`, `GRC_vs_Outgroup.tsv`). 
> Run Step 10 before executing Step 9.

```bash
python 09_alien_index.py
```

AI values are then further processed and visualised using the accompanying R script:

```r
# Run in outputs/alien_index/
source("alien_index.R")
```

| Output | Description |
|--------|-------------|
| `outputs/alien_index/Expressed_Genes_Alien_Index_summary.csv` | Alien Index values for all confidently expressed GRC-linked transcripts |
| `outputs/alien_index/Alien_Index_summary.txt` | Summary of Alien Index classification results |
| `outputs/alien_index/Alien_Index_summary_Supplementry_version.xlsx` | Supplementary-formatted version of the Alien Index summary table |

> ***Note:*** Results are used to produce Figure 3E. All four expressed TE-like transcripts 
> showed greater relative similarity to the *B. coprophila* core proteome than to cecidomyiid 
> or mosquito counterparts, suggesting these loci originated via translocation from the sciarid 
> core genome onto the GRC. Ancestry of insect-like transcripts was generally less resolved, 
> with most showing low similarity to all three proteomes.

---

### Step 10 — Three-way BLASTp: sciarid vs. cecidomyiid vs. mosquito outgroup

To generate the BLASTp inputs required for the Alien Index calculation in Step 9, we perform 
BLASTp searches of expressed GRC protein sequences against three reference proteomes:

- **Sciarid core**: *B. coprophila* core chromosome (non-GRC) proteome
- **Cecidomyiid**: *Aphidoletes aphidimyza* proteome — the closest phylogenetic neighbour 
  to cecidomyiid-derived GRC-linked genes (Hodson et al. 2026)
- **Dipteran outgroup**: *Anopheles gambiae* proteome (GCF_943734735.2)

```bash
10_sciarid_v_cecidomyiid_v_mosquito_BLASTp.sh
```

For each search, the best hit (e-value ≤ 10, `max_target_seqs = 10`) is retained.

| Output | Description |
|--------|-------------|
| `outputs/alien_index/GRC_vs_SciaridCore.tsv` | BLASTp hits of expressed GRC proteins against the *B. coprophila* core proteome |
| `outputs/alien_index/GRC_vs_Cecidomyiid.tsv` | BLASTp hits of expressed GRC proteins against the *A. aphidimyza* proteome |
| `outputs/alien_index/GRC_vs_Outgroup.tsv` | BLASTp hits of expressed GRC proteins against the *A. gambiae* proteome |

> ***Note:*** Transcripts whose outgroup (*A. gambiae*) bitscore exceeds 80% of their best 
> within-clade bitscore are flagged as potentially reflecting broad dipteran conservation 
> rather than lineage-specific ancestry, reducing confidence in ancestral inference. For 
> example, g18444 shows comparable bitscores against both the cecidomyiid and mosquito 
> proteomes, suggesting its apparent cecidomyiid affinity may reflect broad sequence 
> conservation across Diptera rather than a GRC origin through hybridisation.

---

### Step 11 — Validate Alien Index using the complete GRC proteome

To validate the Alien Index as a reliable method for inferring GRC-linked gene ancestry, 
we extend the analysis from the 10 confidently expressed loci to the entire annotated 
*B. coprophila* GRC proteome (22,133 transcripts). This tests whether the approach 
recovers the expected mixed cecidomyiid/sciarid origin of GRC gene content documented 
in Hodson et al. (2026).

```bash
11_complete_GRC_AI.sh
```

Results are processed using:

```r
# Run in outputs/alien_index/
source("complete_GRC_proteome_alien_index.R")
```

| Output | Description |
|--------|-------------|
| `outputs/alien_index/Full_GRC_proteome_AI_summary.csv` | Alien Index values for all 22,133 annotated *B. coprophila* GRC proteome transcripts |

> ***Note:*** Of the 22,133 GRC transcripts, 8,090 (36.6%) were classified as 
> cecidomyiid-like (AI > 10), 7,373 (33.3%) as sciarid-like (AI < −10), and 6,670 (30.1%) 
> as ambiguous. These proportions are broadly consistent with the mixed phylogenetic origin 
> of GRC gene content identified by Hodson et al. (2026) using BUSCO-based ancestry 
> inference, validating the Alien Index approach before application to the smaller set of 
> confidently expressed loci (Supplementary Figure S11).
```
