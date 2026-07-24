library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)
library(dplyr)

box_theme <- theme_bw() + 
  theme(plot.title = element_text(color="black", size=22, face="bold"),
        panel.background = element_rect(fill="white"),
        axis.title.x = element_text(color="black", size=15),
        axis.text.x = element_text(size = 15, hjust = 1),
        axis.title.y = element_text(color="black", size=15),
        axis.text.y = element_text(size = 15), 
        legend.position = "right", 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))

### ---------------------------------------------------------------------------------------###
### PLOT A — All genes expressed per chromosome per stage (genome-wide)
### ---------------------------------------------------------------------------------------###
chromosomes <- c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X", "SUPER_GRC1", "SUPER_GRC2")
GRCs        <- c("SUPER_GRC1", "SUPER_GRC2")

# Load gene counts from GFF3
gtf <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Annotations\\bcop_core_GRC.gff3",
                comment = "#", col_names = FALSE)
colnames(gtf)[1:9] <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Get gene no. per scaffold 
genes_per_scaffold <- gtf %>%
  filter(feature == "gene") %>%
  mutate(gene_id = str_extract(attribute, "ID=[^;]+") %>% str_remove("ID=")) %>%
  distinct(seqname, gene_id) %>%
  filter(seqname %in% chromosomes) %>%
  count(seqname, name = "Num_genes") %>%
  mutate(Chromosome = fct_relevel(seqname, chromosomes))
genes_per_scaffold

# Load TPM + expression status for all genes 
TPM_genes <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\02_intergenic_TPM\\TPM_genes.tsv")
TPM_genes 

# Apply transcription threshold  + reproducibility filter (2/3 expression across libraries per sex) across all chromosomes
# Use same logic as GRC filter but genome-wide, per sex/stage/tissue
expected_libraries_all <- TPM_genes %>%
  filter(Chromosome %in% chromosomes) %>%
  distinct(Sex, Stage, Tissue, Sample) %>%
  group_by(Sex, Stage, Tissue) %>%
  summarise(n_libraries = n(), .groups = "drop") %>%
  mutate(n_required = ceiling(2/3 * n_libraries))

gene_counts_all <- TPM_genes %>%
  filter(
    Chromosome %in% chromosomes,
    Expression_Status_BH == "Active"
  ) %>%
  group_by(Gene, Chromosome, Sex, Stage, Tissue) %>%
  summarise(n_above_threshold = n(), .groups = "drop")

# Gene only passes if n above expressed TPM threshold > N-required (2/3 x n)
passing_combinations_all <- gene_counts_all %>%
  left_join(expected_libraries_all, by = c("Sex", "Stage", "Tissue")) %>%
  mutate(passes = n_above_threshold >= n_required) %>%
  filter(passes) %>%
  distinct(Gene, Chromosome, Stage)

# Quick count of core genes v GRC gene expressed
core_genes_expressd <- passing_combinations_all %>%
  filter(Chromosome %in% c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X")) %>%
  distinct(Gene) %>%
  nrow()
cat("Core genes expressed:", core_genes_expressd)

GRC_genes_expressd <- passing_combinations_all %>%
  filter(Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2")) %>%
  distinct(Gene) %>%
  nrow()
cat("GRC genes expressed:", GRC_genes_expressd)

# No of time GRC genes are expressed
passing_combinations_all %>%
  filter(Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2")) %>%
  nrow()

# Breakdown per chromosome if needed
passing_combinations_all %>%
  filter(Chromosome %in% c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X", "SUPER_GRC1", "SUPER_GRC2")) %>%
  group_by(Chromosome) %>%
  summarise(n_expressed_genes = n_distinct(Gene), .groups = "drop")

# Count expressed genes per chromosome per stage
counts_per_sample <- passing_combinations_all %>%
  group_by(Chromosome, Stage) %>%
  summarise(Num_genes = n_distinct(Gene), .groups = "drop") %>%
  mutate(Chromosome = fct_relevel(Chromosome, chromosomes))
counts_per_sample

# Pivot and join total gene counts
wide_expr <- counts_per_sample %>%
  pivot_wider(names_from = Stage, values_from = Num_genes, values_fill = 0)

combined_counts <- genes_per_scaffold %>%
  rename(Total = Num_genes) %>%
  left_join(wide_expr, by = "Chromosome")
combined_counts

combined_long <- combined_counts %>%
  select(Chromosome, Total, `0-4h`, `4-8h`, `late-larva-early-pupa`, adult) %>%
  pivot_longer(cols = -Chromosome, names_to = "Stage", values_to = "Count") %>%
  mutate(
    Stage      = factor(Stage, levels = c("Total", "0-4h", "4-8h", "late-larva-early-pupa", "adult")),
    Chromosome = fct_relevel(Chromosome, chromosomes)
  )

plot_compare <- ggplot(combined_long, aes(x = Chromosome, y = Count, fill = Stage)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(x = "Chromosome", y = "Gene Count", title = "A") +
  box_theme +
  scale_fill_manual(
    values = c(
      "Total"                 = "grey70",
      "0-4h"                  = "#f3c6f4",
      "4-8h"                  = "#e088ec",
      "late-larva-early-pupa" = "#c043d6",
      "adult"                 = "#7a0177"
    ),
    labels = c(
      "Total"                 = "Total number of genes",
      "0-4h"                  = "Genes expressed (0–4h embryo)",
      "4-8h"                  = "Genes expressed (4–8h embryo)",
      "late-larva-early-pupa" = "Genes expressed (larval/pupa)",
      "adult"                 = "Genes expressed (adult)"
    ),
    name = NULL
  ) +
  scale_x_discrete(labels = c(
    "SUPER_1" = "I", "SUPER_2" = "II", "SUPER_3" = "III",
    "SUPER_X" = "X", "SUPER_GRC1" = "GRC1", "SUPER_GRC2" = "GRC2"
  ))

plot_compare

### ---------------------------------------------------------------------------------------###
### PLOT B — GRC expressed genes per stage (all passing reproducibility, before extra filters)
### ---------------------------------------------------------------------------------------###

GRC_counts_long <- passing_combinations_all %>%
  filter(Chromosome %in% GRCs) %>%
  group_by(Chromosome, Stage) %>%
  summarise(Num_genes = n_distinct(Gene), .groups = "drop") %>%
  mutate(
    Stage      = factor(Stage, levels = c("0-4h", "4-8h", "late-larva-early-pupa", "adult")),
    Chromosome = fct_relevel(Chromosome, GRCs)
  )

plot_GRC_flat <- ggplot(GRC_counts_long, aes(x = Chromosome, y = Num_genes, fill = Stage)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(x = "Chromosome", y = "Number of Expressed Genes", title = "B") +
  box_theme +
  scale_fill_manual(
    values = c(
      "0-4h"                  = "#f3c6f4",
      "4-8h"                  = "#e088ec",
      "late-larva-early-pupa" = "#c043d6",
      "adult"                 = "#7a0177"
    ),
    labels = c(
      "0-4h"                  = "Genes expressed (0–4h embryo)",
      "4-8h"                  = "Genes expressed (4–8h embryo)",
      "late-larva-early-pupa" = "Genes expressed (larval/pupa)",
      "adult"                 = "Genes expressed (adult)"
    ),
    name = NULL
  ) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("SUPER_GRC1" = "GRC1", "SUPER_GRC2" = "GRC2"))

plot_GRC_flat

### ---------------------------------------------------------------------------------------###
### PLOT C — Final filtered GRC genes (germline-only, low alignment) per stage
### ---------------------------------------------------------------------------------------###

Final_expressed_GRC_genes <- Final_expressed_GRC_genes %>%
  filter(gene_id != c("g19121", "g19161")) # Remove bacterial contaminants 

GRC_filtered_counts_long <- Final_expressed_GRC_genes %>%
  group_by(scaffold, development_stage) %>%
  summarise(Num_genes = n_distinct(gene_id), .groups = "drop") %>%
  rename(Chromosome = scaffold, Stage = development_stage) %>%
  mutate(
    Stage      = factor(Stage, levels = c("0-4h", "4-8h", "late-larva-early-pupa", "adult")),
    Chromosome = fct_relevel(Chromosome, GRCs)
  )


plot_GRC_corrected <- ggplot(GRC_filtered_counts_long, aes(x = Chromosome, y = Num_genes, fill = Stage)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(x = "Chromosome", y = "Number of Expressed Genes", title = "C") +
  box_theme +
  scale_fill_manual(
    values = c(
      "0-4h"                  = "#f3c6f4",
      "4-8h"                  = "#e088ec",
      "late-larva-early-pupa" = "#c043d6",
      "adult"                 = "#7a0177"
    ),
    labels = c(
      "0-4h"                  = "Genes expressed (0–4h embryo)",
      "4-8h"                  = "Genes expressed (4–8h embryo)",
      "late-larva-early-pupa" = "Genes expressed (larval/pupa)",
      "adult"                 = "Genes expressed (adult)"
    ),
    name = NULL
  ) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("SUPER_GRC1" = "GRC1", "SUPER_GRC2" = "GRC2"))

plot_GRC_corrected

### ---------------------------------------------------------------------------------------###
### COMBINE
### ---------------------------------------------------------------------------------------###

# Get shared y max for B and C
y_max_BC <- max(
  max(GRC_counts_long$Num_genes),
  max(GRC_filtered_counts_long$Num_genes)
) + 2

p_final <- plot_compare / (
  (plot_GRC_flat      + scale_y_continuous(limits = c(0, y_max_BC))) +
    (plot_GRC_corrected + scale_y_continuous(limits = c(0, y_max_BC)))
)

p_final


