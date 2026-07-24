library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

### ---------------------------------------------------------------------------------------###
### Build df_tpms and df_summary_tpms from Expressed_GRC_genes_summary
### ---------------------------------------------------------------------------------------###
setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\03_expressed_GRC_genes\\outputs")
home <- getwd()

Expressed_GRC_genes_summary <- read_tsv("Expressed_GRC_genes_summary.tsv")
Expressed_GRC_genes_summary
df_tpms <- Expressed_GRC_genes_summary %>%
  rename(
    gene_id           = gene_id,
    development_stage = development_stage
  ) %>%
  separate(`mean_germ_TPM/mean_soma_TPM`,
           into = c("mean_germ_TPM", "mean_soma_TPM"),
           sep = "_", convert = TRUE) %>%
  mutate(
    development_stage = factor(development_stage,
                               levels = c("0-4h", "4-8h", "late-larva-early-pupa", "adult"))
  )
df_tpms
df_summary_tpms <- df_tpms %>%
  group_by(gene_id, development_stage) %>%
  summarise(
    mean_germ_TPM   = max(mean_germ_TPM, na.rm = TRUE),
    mean_soma_TPM   = max(mean_soma_TPM, na.rm = TRUE),
    alignment_score = max((Coverage * `%Identity`) / 100, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(alignment_score = ifelse(is.infinite(alignment_score) | is.na(alignment_score), 0, alignment_score))

# ---- Build heatmap long-format ----
# Pre-GRC stages: average germ/soma into one column
# Late-larva/adult: keep germ and soma separate
df_heat <- df_summary_tpms %>%
  pivot_longer(cols = c(mean_germ_TPM, mean_soma_TPM),
               names_to = "library_type", values_to = "TPM") %>%
  mutate(
    col_id = case_when(
      development_stage %in% c("0-4h", "4-8h") ~ as.character(development_stage),
      library_type == "mean_germ_TPM"            ~ paste0(development_stage, "\nGermline"),
      library_type == "mean_soma_TPM"             ~ paste0(development_stage, "\nSomatic")
    )
  ) %>%
  group_by(gene_id, col_id) %>%
  summarise(log_TPM = log1p(mean(TPM, na.rm = TRUE)), .groups = "drop")
df_heat
# ---- Build matrix & hierarchically cluster rows + columns ----
heat_mat <- df_heat %>%
  pivot_wider(names_from = col_id, values_from = log_TPM, values_fill = 0) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()
heat_mat

# Enforce fixed column order (no column clustering)
logical_col_order <- c(
  "0-4h",
  "4-8h",
  "late-larva-early-pupa\nGermline",
  "late-larva-early-pupa\nSomatic",
  "adult\nGermline",
  "adult\nSomatic"
)

heat_mat <- heat_mat[, logical_col_order[logical_col_order %in% colnames(heat_mat)]]

# Cluster rows only, keep columns fixed
row_order <- rownames(heat_mat)[hclust(dist(heat_mat))$order]
col_order <- logical_col_order[logical_col_order %in% colnames(heat_mat)]  

# ---- Apply cluster ordering to plot data ----
df_heat_plot <- df_heat %>%
  mutate(
    gene_id = factor(gene_id, levels = row_order),
    col_id  = factor(col_id,  levels = col_order)
  )

# ---- Stage group annotation (top bar above heatmap) ----
stage_group_map <- c(
  "0-4h"                              = "pre-GRC elimination",
  "4-8h"                              = "pre-GRC elimination",
  "late-larva-early-pupa\nGermline"   = "Germline",
  "late-larva-early-pupa\nSomatic"    = "Somatic",
  "adult\nGermline"                   = "Germline",
  "adult\nSomatic"                    = "Somatic"
)
stage_colors <- c(
  "pre-GRC elimination"      = "#AED6F1",
  "Germline"  = "#A9DFBF",
  "Somatic"                    = "#ffccaa"
)

# ---- Column annotation bar (top) ----
df_col_ann <- data.frame(col_id = factor(col_order, levels = col_order)) %>%
  mutate(stage_group = factor(stage_group_map[as.character(col_id)],
                              levels = names(stage_colors)))

p_col_ann <- ggplot(df_col_ann, aes(x = col_id, y = 1, fill = stage_group)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(values = stage_colors, name = "Stage") +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none"   # legend shown in main plot instead
  )

# ---- Main heatmap ----
p_heat <- ggplot(df_heat_plot, aes(x = col_id, y = gene_id, fill = log_TPM)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#2C3E90",
                      name = "log(TPM + 1)") +
  scale_x_discrete(position = "bottom") +  # <-- moved to bottom
  labs(x = NULL, y = "Gene ID") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),  # hjust=1 for bottom angled labels
    axis.text.y  = element_text(size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  )

# ---- Sex annotation strip (right side) ----
df_sex <- df_tpms %>%
  group_by(gene_id) %>%
  summarise(
    all_sexes = list(unique(sex)),
    .groups = "drop"
  ) %>%
  mutate(
    sex = case_when(
      sapply(all_sexes, function(x) all(x == "male"))        ~ "male",
      sapply(all_sexes, function(x) all(x == "female"))      ~ "female",
      TRUE                                                    ~ "both-sexes"
    )
  ) %>%
  select(gene_id, sex) %>%
  mutate(gene_id = factor(gene_id, levels = row_order))

df_sex

sex_colors <- c(
  "male"       = "#4A90D9",   # Steel blue
  "female"     = "#E8527A",   # Rose pink
  "both-sexes" = "#8B6DB5"    # Muted purple (blend of the two)
)

p_sex <- ggplot(df_sex, aes(x = 1, y = gene_id, fill = sex)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = sex_colors, name = "Sex\nExpression") +
  scale_x_continuous(breaks = 1, labels = "Sex") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 10, angle = 45, hjust = 0),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  )
# ---- Alignment score strip (right side) ----
df_align <- df_summary_tpms %>%
  group_by(gene_id) %>%
  summarise(alignment_score = max(alignment_score, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    alignment_score = ifelse(is.infinite(alignment_score), 0, alignment_score),
    gene_id = factor(gene_id, levels = row_order)
  )

p_align <- ggplot(df_align, aes(x = 1, y = gene_id, fill = alignment_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma",
                       name = "Alignment Score\n(%Cov × %Id / 100)",
                       na.value = "white") +
  scale_x_continuous(breaks = 1, labels = "Alignment") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 10, angle = 45, hjust = 0),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  )

# ---- Combine with patchwork ----
# Layout: main heatmap wide, two narrow annotation strips on right
p_final <- p_heat + p_align +
  plot_layout(widths = c(6, 0.4, 0.4), guides = "collect") +
  plot_annotation(
    title    = "GRC-linked gene expression across developmental stages",
    subtitle = "Genes and columns ordered by hierarchical clustering",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  ) &
  theme(legend.position = "bottom")

p_final

# ---- Blank spacer to align top bar with annotation strips ----
p_blank <- ggplot() + theme_void()

# ---- Combine with patchwork ----
# Top row: col annotation bar + two blank spacers (aligning with sex/align strips)
# Bottom row: main heatmap + sex strip + align strip
# ---- Combine: top bar sits directly on heatmap ----
top_row    <- p_col_ann + p_blank + p_blank +
  plot_layout(widths = c(6, 0.4, 0.4))

bottom_row <- p_heat + p_align +
  plot_layout(widths = c(6, 0.4, 0.4), guides = "collect")

p_final <- top_row / bottom_row +
  plot_layout(heights = c(0.06, 1)) +  # slightly taller bar so colour is visible
  plot_annotation(
    title    = "GRC-linked gene expression across developmental stages",
    subtitle = "Genes ordered by hierarchical clustering",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  ) &
  theme(
    legend.position    = "bottom",
    legend.title.position = "top",   # <-- all legend titles above their keys
    legend.title       = element_text(hjust = 0.5)
  )
p_final


# ---- Summarise stage group mid-points for text labels ----
df_col_ann_labels <- df_col_ann %>%
  mutate(x_num = as.integer(col_id)) %>%
  group_by(stage_group) %>%
  summarise(x_mid = mean(x_num), .groups = "drop")

# ---- Combine (no change needed here now) ----
top_row    <- p_col_ann + p_blank + p_blank +
  plot_layout(widths = c(6, 0.4, 0.4))

bottom_row <- p_heat + p_sex + p_align +
  plot_layout(widths = c(6, 0.4, 0.4), guides = "collect")

p_final <- top_row / bottom_row +
  plot_layout(heights = c(0.06, 1)) +
  plot_annotation(
    title    = "GRC-linked gene expression across developmental stages",
    subtitle = "Genes ordered by hierarchical clustering",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  ) &
  theme(
    legend.position       = "bottom",
    legend.title.position = "top",
    legend.title          = element_text(hjust = 0.5)
  )

p_final

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions")
home <- getwd()
#ggsave("GRC_heatmap_unedited.svg", p_final, width = 14, height = max(8, nrow(heat_mat) * 0.25 + 3))
