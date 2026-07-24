# ======================================================================================
# GRC Gene Expression Heatmap — Mean TPM per sex/stage/tissue
# ======================================================================================
# This script visualises expression of actively expressed GRC genes as a heatmap,
# where each column represents the mean log2(TPM) across biological replicates
# within a sex/stage/tissue group. Genes are ordered by hierarchical clustering
# of their binary activity pattern (expressed/not expressed per group), which
# groups genes by where they are expressed rather than by TPM magnitude.
#
# Adult females: gynogenic (B10-12, B16-18) and androgenic (B19-21, B25-27)
# replicates are averaged together into a single "Female Adult" group per tissue,
# giving 6 replicates per adult female tissue group.
#
# Input:  Expressed_GRC_genes — filtered TPM data from GRC_expression_analysis.R
#         blast_df             — BLAST alignment results loaded in same session
# ======================================================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

# ======================================================================================
# STEP 1: Map each sample to its sex/stage/tissue group
# ======================================================================================
# Adult female gynogenic (G) and androgenic (A) replicates are merged into a
# single "Female Adult Germ" and "Female Adult Soma" group respectively,
# averaging across all 6 replicates per tissue.



sample_group_map <- c(
  ME1    = "Male 0-4h",        ME2    = "Male 0-4h",        ME3    = "Male 0-4h",
  FE1    = "Female 0-4h",      FE2    = "Female 0-4h",      FE3    = "Female 0-4h",
  ML1    = "Male 4-8h",        ML2    = "Male 4-8h",        ML3    = "Male 4-8h",
  FL1    = "Female 4-8h",      FL2    = "Female 4-8h",      FL3    = "Female 4-8h",
  Mgerm1 = "Male LL/EP Germ",  Mgerm2 = "Male LL/EP Germ",  Mgerm3 = "Male LL/EP Germ",
  Mbody1 = "Male LL/EP Soma",  Mbody2 = "Male LL/EP Soma",  Mbody3 = "Male LL/EP Soma",
  Fgerm1 = "Female LL/EP Germ",Fgerm2 = "Female LL/EP Germ",Fgerm3 = "Female LL/EP Germ",
  Fbody1 = "Female LL/EP Soma",Fbody2 = "Female LL/EP Soma",Fbody3 = "Female LL/EP Soma",
  B1     = "Male Adult Germ",  B2     = "Male Adult Germ",  B3     = "Male Adult Germ",
  B7     = "Male Adult Soma",  B8     = "Male Adult Soma",  B9     = "Male Adult Soma",
  # Gynogenic females — merged with androgenic into single Female Adult group
  B10    = "Female Adult Germ", B11   = "Female Adult Germ", B12   = "Female Adult Germ",
  B16    = "Female Adult Soma", B17   = "Female Adult Soma", B18   = "Female Adult Soma",
  # Androgenic females — merged with gynogenic into single Female Adult group
  B19    = "Female Adult Germ", B20   = "Female Adult Germ", B21   = "Female Adult Germ",
  B25    = "Female Adult Soma", B26   = "Female Adult Soma", B27   = "Female Adult Soma"
)

# Fixed column order — embryo stages first, then larval, then adult;
# male before female; germ before soma
group_col_order <- c(
  "Male 0-4h",
  "Female 0-4h",
  "Male 4-8h",
  "Female 4-8h",
  "Male LL/EP Germ",
  "Female LL/EP Germ",
  "Male LL/EP Soma",
  "Female LL/EP Soma",
  "Male Adult Germ",
  "Female Adult Germ",
  "Male Adult Soma",
  "Female Adult Soma"
)

# ======================================================================================
# STEP 2: Compute mean log2(TPM) per gene per group
# ======================================================================================
# Average log2(TPM) across all replicates within each group.
# For adult females this averages across 6 replicates (3 gynogenic + 3 androgenic).

heat_df_grouped <- Expressed_GRC_genes %>%
  select(Gene, Sample, log2_TPM) %>%
  mutate(group = sample_group_map[Sample]) %>%
  group_by(Gene, group) %>%
  summarise(mean_log2_TPM = mean(log2_TPM, na.rm = TRUE), .groups = "drop")

# ======================================================================================
# STEP 3: Build wide matrix and cluster rows by activity pattern
# ======================================================================================

heat_mat <- heat_df_grouped %>%
  pivot_wider(names_from = group, values_from = mean_log2_TPM, values_fill = 0) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

# Keep only groups present in the data, in the defined biological order
col_order <- group_col_order[group_col_order %in% colnames(heat_mat)]
heat_mat  <- heat_mat[, col_order]

# Cluster rows by binary activity pattern (expressed/not expressed per group)
# rather than by TPM magnitude — this groups genes by where they are expressed
# rather than how highly, making expression patterns easier to interpret visually
activity_mat <- heat_mat
activity_mat[activity_mat > 0] <- 1
row_order <- rownames(activity_mat)[hclust(dist(activity_mat), method = "ward.D2")$order]

# ======================================================================================
# STEP 4: Build long-format plot data
# ======================================================================================

df_heat_plot <- heat_df_grouped %>%
  mutate(
    Gene  = factor(Gene,  levels = row_order),
    group = factor(group, levels = col_order)
  )

# ======================================================================================
# STEP 5: Stage annotation bar (top of heatmap)
# ======================================================================================
# Columns are coloured by developmental context:
#   pre-GRC elimination — embryonic stages before germline genome elimination
#   Germline            — germline tissue at larval/adult stages
#   Somatic             — somatic tissue at larval/adult stages

stage_group_map <- c(
  "Male 0-4h"         = "pre-GRC elimination",
  "Female 0-4h"       = "pre-GRC elimination",
  "Male 4-8h"         = "pre-GRC elimination",
  "Female 4-8h"       = "pre-GRC elimination",
  "Male LL/EP Germ"   = "Germline",
  "Female LL/EP Germ" = "Germline",
  "Male LL/EP Soma"   = "Somatic",
  "Female LL/EP Soma" = "Somatic",
  "Male Adult Germ"   = "Germline",
  "Female Adult Germ" = "Germline",
  "Male Adult Soma"   = "Somatic",
  "Female Adult Soma" = "Somatic"
)

stage_colors <- c(
  "pre-GRC elimination" = "#AED6F1",
  "Germline"            = "#A9DFBF",
  "Somatic"             = "#ffccaa"
)

df_col_ann <- data.frame(
  group = factor(col_order, levels = col_order)
) %>%
  mutate(stage_group = factor(stage_group_map[as.character(group)],
                              levels = names(stage_colors)))

p_col_ann <- ggplot(df_col_ann, aes(x = group, y = 1, fill = stage_group)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(values = stage_colors, name = "Stage") +
  theme_void() +
  theme(legend.position = "none")

# ======================================================================================
# STEP 6: Main heatmap
# ======================================================================================
p_heat <- ggplot(df_heat_plot, aes(x = group, y = Gene, fill = mean_log2_TPM)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#EDF0FA", high = "#2C3E90",
                      name = "mean log2(TPM)") +
  labs(x = NULL, y = "Gene ID") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y     = element_text(size = 7),
    panel.grid      = element_blank(),
    legend.position = "bottom"
  )

# ======================================================================================
# STEP 7: Alignment score annotation strip (right side)
# ======================================================================================
# Alignment score = (Coverage x %Identity) / 100, derived from BLAST search of
# GRC genes against the core genome. A score of 0 indicates no significant hit;
# higher scores suggest potential homology to core-genome loci.
blast_df <- read_tsv("03_expressed_GRC_genes\\GRC_BLAST_table.tsv") %>%
  mutate(
    `%Identity` = replace_na(`%Identity`, 0),
    Coverage    = replace_na(Coverage, 0)
  )
blast_df

alignment_scores <- data.frame(gene_id = row_order) %>%
  left_join(
    blast_df %>%
      mutate(alignment_score = (Coverage * `%Identity`) / 100) %>%
      group_by(gene_id) %>%
      summarise(alignment_score = max(alignment_score, na.rm = TRUE), .groups = "drop") %>%
      mutate(alignment_score = ifelse(is.infinite(alignment_score) | is.na(alignment_score),
                                      0, alignment_score)),
    by = "gene_id"
  ) %>%
  mutate(
    alignment_score = replace_na(alignment_score, 0),
    gene_id = factor(gene_id, levels = row_order)
  )
alignment_scores

# Genes in the heatmap but absent from blast_df
missing_from_blast <- setdiff(row_order, blast_df$gene_id)
missing_from_blast

p_align <- ggplot(alignment_scores, aes(x = 1, y = gene_id, fill = alignment_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option   = "magma",
                       name     = "Alignment Score\n(%Cov × %Id / 100)",
                       na.value = "white") +
  scale_x_continuous(breaks = 1, labels = "Alignment") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y     = element_blank(),
    axis.text.x     = element_text(size = 9, angle = 45, hjust = 0),
    panel.grid      = element_blank(),
    legend.position = "bottom"
  )
p_align
# ======================================================================================
# STEP 8: Combine with patchwork
# ======================================================================================

p_blank <- ggplot() + theme_void()

top_row    <- p_col_ann + p_blank +
  plot_layout(widths = c(10, 0.4))

bottom_row <- p_heat + p_align +
  plot_layout(widths = c(10, 0.4), guides = "collect")

p_final_grouped <- top_row / bottom_row +
  plot_layout(heights = c(0.05, 1)) +
  plot_annotation(
    title    = "Expressed GRC genes — mean TPM per sex/stage/tissue group",
    subtitle = "Genes ordered by activity pattern clustering | Adult females averaged across gynogenic and androgenic",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  ) &
  theme(
    legend.position       = "bottom",
    legend.title.position = "top",
    legend.title          = element_text(hjust = 0.5)
  )

p_final_grouped

# ======================================================================================
# STEP 9: Save
# ======================================================================================

