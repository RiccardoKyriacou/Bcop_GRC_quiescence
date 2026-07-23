# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)

################
# Expressed GRC gene analysis
################

TPM_genes <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\02_intergenic_TPM\\TPM_genes.tsv")
TPM_genes 

Expressed_GRC_genes <- TPM_genes %>% 
  filter(Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2"),
         Expression_Status_BH == "Active")
Expressed_GRC_genes

# Step 1: Define expected n per sex/stage/tissue from the full design
expected_libraries <- TPM_genes %>%
  filter(Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2")) %>%
  distinct(Sex, Stage, Tissue, Sample) %>%
  group_by(Sex, Stage, Tissue) %>%
  summarise(n_libraries = n(), .groups = "drop") %>%
  mutate(n_required = ceiling(2/3 * n_libraries))

expected_libraries

# Step 2: Count how many libraries each gene passes threshold in per sex/stage/tissue
gene_counts <- TPM_genes %>%
  filter(
    Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2"),
    Expression_Status_BH == "Active"
  ) %>%
  group_by(Gene, Sex, Stage, Tissue) %>%
  summarise(n_above_threshold = n(), .groups = "drop")
gene_counts

# Step 3: Join expected counts and apply threshold
reproducibility_filter <- gene_counts %>%
  left_join(expected_libraries, by = c("Sex", "Stage", "Tissue")) %>%
  mutate(passes = n_above_threshold >= n_required)
reproducibility_filter

# Step 4: Keep gene/stage/tissue combinations that pass — not just the gene
passing_combinations <- reproducibility_filter %>%
  filter(passes) %>%
  select(Gene, Sex, Stage, Tissue)

# Step 5: Final filter — only keep rows that belong to a passing combination
Expressed_GRC_genes <- TPM_genes %>%
  filter(
    Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2"),
    Expression_Status_BH == "Active"
  ) %>%
  inner_join(passing_combinations, by = c("Gene", "Sex", "Stage", "Tissue"))

Expressed_GRC_genes
### ---------------------------------------------------------------------------------------###
### Build Expressed_GRC_genes summary table
### ---------------------------------------------------------------------------------------###

# Classify each tissue as germ or soma
tissue_type_map <- c(
  "germ_embryo"  = "germ",
  "germ"         = "germ",
  "germ_testes"  = "germ",
  "germ_ovaries" = "germ",
  "soma_carcass" = "soma"
)

# Build summary per gene x stage x sex
expressed_summary <- Expressed_GRC_genes %>%
  mutate(tissue_type = tissue_type_map[Tissue]) %>%
  group_by(Gene, Chromosome, Stage, Sex, tissue_type) %>%
  summarise(
    n_libraries  = n(),
    mean_TPM     = mean(TPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Pivot to get germ and soma side by side
  pivot_wider(
    names_from  = tissue_type,
    values_from = c(n_libraries, mean_TPM),
    values_fill = 0
  ) %>%
  # Rename for clarity
  rename_with(~ gsub("n_libraries_", "n_", .), starts_with("n_libraries_")) %>%
  rename_with(~ gsub("mean_TPM_", "mean_TPM_", .), starts_with("mean_TPM_"))
expressed_summary

# Summarise sex per gene x stage (male/female/both-sexes)
sex_summary <- expressed_summary %>%
  group_by(Gene, Chromosome, Stage) %>%
  summarise(
    sex = case_when(
      all(Sex == "male")   ~ "male",
      all(Sex == "female") ~ "female",
      TRUE                 ~ "both-sexes"
    ),
    n_germ       = sum(n_germ,       na.rm = TRUE),
    n_soma       = sum(n_soma,       na.rm = TRUE),
    mean_TPM_germ = mean(mean_TPM_germ, na.rm = TRUE),
    mean_TPM_soma = mean(mean_TPM_soma, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # tissue category: germ only, soma only, or both
    tissue_category = case_when(
      n_germ > 0 & n_soma > 0 ~ "both",
      n_germ > 0              ~ "germ",
      TRUE                    ~ "soma"
    ),
    `germ/soma`               = paste0(n_germ, "_", n_soma),
    `mean_germ_TPM/mean_soma_TPM` = paste0(round(mean_TPM_germ, 2), "_", round(mean_TPM_soma, 2))
  ) %>%
  select(
    gene_id   = Gene,
    scaffold  = Chromosome,
    tissue_category,
    development_stage = Stage,
    sex,
    `germ/soma`,
    `mean_germ_TPM/mean_soma_TPM`
  )

# Join BLAST alignment data
blast_df <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\03_expressed_GRC_genes\\GRC_BLAST_table.tsv") %>%
  mutate(
    `%Identity` = replace_na(`%Identity`, 0),
    Coverage    = replace_na(Coverage, 0)
  )

blast_summary <- blast_df %>%
  mutate(`%Identity` = replace_na(`%Identity`, 0),
         Coverage    = replace_na(Coverage, 0)) %>%
  group_by(gene_id) %>%
  slice_max(order_by = (Coverage * `%Identity`) / 100, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene_id, Hit_id, Coverage, `%Identity`, `E-Value`, Hit_length, Mismatches)

# Final join
Expressed_GRC_genes_summary <- sex_summary %>%
  left_join(blast_summary, by = "gene_id") %>%
  arrange(gene_id, development_stage)

Expressed_GRC_genes_summary

# Final filter
GRC_germline_only_no_hit <- Expressed_GRC_genes_summary %>%
  # Split the TPM column to get soma TPM
  separate(`mean_germ_TPM/mean_soma_TPM`, 
           into = c("mean_germ_TPM", "mean_soma_TPM"), 
           sep = "_", convert = TRUE) %>%
  group_by(gene_id) %>%
  filter(
    # Soma TPM is 0 across all stages
    all(mean_soma_TPM == 0),
    # Alignment score < 70 across all stages
    all((Coverage * `%Identity`) / 100 < 70)
  ) %>%
  ungroup() %>%
  distinct(gene_id)

GRC_germline_only_no_hit

Final_expressed_GRC_genes <- Expressed_GRC_genes_summary %>%
  filter(gene_id %in% GRC_germline_only_no_hit$gene_id)
Final_expressed_GRC_genes

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\03_expressed_GRC_genes\\outputs")
home <- getwd()
write.table(Expressed_GRC_genes,
            file = file.path(home, "Expressed_GRC_genes.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
write.table(Expressed_GRC_genes_summary,
            file = file.path(home, "Expressed_GRC_genes_summary.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
write.table(Final_expressed_GRC_genes,
            file = file.path(home, "Corrected_expressed_GRC_genes.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


