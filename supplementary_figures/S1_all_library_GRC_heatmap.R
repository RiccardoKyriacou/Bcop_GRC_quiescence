

### ---------------------------------------------------------------------------------------###
### HEATMAP: Expressed GRC genes across all replicates
### ---------------------------------------------------------------------------------------###

# Step 1: Define sample column order (logical grouping by sex/stage/tissue)
sample_order <- c(
  # Embryo male
  "ME1", "ME2", "ME3",
  # Embryo female
  "FE1", "FE2", "FE3",
  # 4-8h male
  "ML1", "ML2", "ML3",
  # 4-8h female
  "FL1", "FL2", "FL3",
  # Late larva male germ
  "Mgerm1", "Mgerm2", "Mgerm3",
  # Late larva male soma
  "Mbody1", "Mbody2", "Mbody3",
  # Late larva female germ
  "Fgerm1", "Fgerm2", "Fgerm3",
  # Late larva female soma
  "Fbody1", "Fbody2", "Fbody3",
  # Adult male germ
  "B1", "B2", "B3",
  # Adult male soma
  "B7", "B8", "B9",
  # Adult female germ (gynogenic)
  "B10", "B11", "B12",
  # Adult female soma (gynogenic)
  "B16", "B17", "B18",
  # Adult female germ (androgenic)
  "B19", "B20", "B21",
  # Adult female soma (androgenic)
  "B25", "B26", "B27"
)

# Step 2: Build wide matrix — one row per gene, one column per sample
# Step 2: Build wide matrix — average duplicates first
heat_df <- Expressed_GRC_genes %>%
  select(Gene, Sample, log2_TPM) %>%
  group_by(Gene, Sample) %>%
  summarise(log2_TPM = mean(log2_TPM, na.rm = TRUE), .groups = "drop") %>%  # handles duplicates
  pivot_wider(names_from = Sample, values_from = log2_TPM, values_fill = 0) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

# Keep only samples present in data, in logical order
col_order <- sample_order[sample_order %in% colnames(heat_df)]
heat_df   <- heat_df[, col_order]

# Step 3: Cluster rows (genes) only
row_order <- rownames(heat_df)[hclust(dist(heat_df))$order]
# Keep only samples present in data, in logical order
col_order <- sample_order[sample_order %in% colnames(heat_df)]
heat_df   <- heat_df[, col_order]

# Step 3: Cluster rows (genes) only
row_order <- rownames(heat_df)[hclust(dist(heat_df))$order]

# Step 4: Build long format for ggplot
df_heat_plot <- Expressed_GRC_genes %>%
  select(Gene, Sample, log2_TPM) %>%
  mutate(
    Gene   = factor(Gene,   levels = row_order),
    Sample = factor(Sample, levels = col_order)
  )

### ---------------------------------------------------------------------------------------###
### COLUMN ANNOTATION BAR
### ---------------------------------------------------------------------------------------###

# Map each sample to its group for the top annotation bar
sample_group_map <- c(
  ME1="Male 0-4h", ME2="Male 0-4h", ME3="Male 0-4h",
  FE1="Female 0-4h", FE2="Female 0-4h", FE3="Female 0-4h",
  ML1="Male 4-8h", ML2="Male 4-8h", ML3="Male 4-8h",
  FL1="Female 4-8h", FL2="Female 4-8h", FL3="Female 4-8h",
  Mgerm1="Male LL/EP Germ", Mgerm2="Male LL/EP Germ", Mgerm3="Male LL/EP Germ",
  Mbody1="Male LL/EP Soma", Mbody2="Male LL/EP Soma", Mbody3="Male LL/EP Soma",
  Fgerm1="Female LL/EP Germ", Fgerm2="Female LL/EP Germ", Fgerm3="Female LL/EP Germ",
  Fbody1="Female LL/EP Soma", Fbody2="Female LL/EP Soma", Fbody3="Female LL/EP Soma",
  B1="Male Adult Germ", B2="Male Adult Germ", B3="Male Adult Germ",
  B7="Male Adult Soma", B8="Male Adult Soma", B9="Male Adult Soma",
  B10="Female Adult Germ (G)", B11="Female Adult Germ (G)", B12="Female Adult Germ (G)",
  B16="Female Adult Soma (G)", B17="Female Adult Soma (G)", B18="Female Adult Soma (G)",
  B19="Female Adult Germ (A)", B20="Female Adult Germ (A)", B21="Female Adult Germ (A)",
  B25="Female Adult Soma (A)", B26="Female Adult Soma (A)", B27="Female Adult Soma (A)"
)

group_colors <- c(
  "Male 0-4h"            = "#AED6F1",
  "Female 0-4h"          = "#F1948A",
  "Male 4-8h"            = "#85C1E9",
  "Female 4-8h"          = "#EC7063",
  "Male LL/EP Germ"      = "#5DADE2",
  "Male LL/EP Soma"      = "#2E86C1",
  "Female LL/EP Germ"    = "#E74C3C",
  "Female LL/EP Soma"    = "#922B21",
  "Male Adult Germ"      = "#1A5276",
  "Male Adult Soma"      = "#154360",
  "Female Adult Germ (G)"= "#C0392B",
  "Female Adult Soma (G)"= "#7B241C",
  "Female Adult Germ (A)"= "#F1948A",
  "Female Adult Soma (A)"= "#CD6155"
)

df_col_ann <- data.frame(
  Sample = factor(col_order, levels = col_order)
) %>%
  mutate(group = factor(sample_group_map[as.character(Sample)],
                        levels = names(group_colors)))

p_col_ann <- ggplot(df_col_ann, aes(x = Sample, y = 1, fill = group)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(values = group_colors, name = "Sample Group") +
  theme_void() +
  theme(legend.position = "none")

### ---------------------------------------------------------------------------------------###
### MAIN HEATMAP
### ---------------------------------------------------------------------------------------###

p_heat <- ggplot(df_heat_plot, aes(x = Sample, y = Gene, fill = log2_TPM)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "#2C3E90", name = "log2(TPM + 1)") +
  labs(x = NULL, y = "Gene ID") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y  = element_text(size = 7),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  )

### ---------------------------------------------------------------------------------------###
### ALIGNMENT SCORE STRIP
### ---------------------------------------------------------------------------------------##

blast_df <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\03_expressed_GRC_genes\\GRC_BLAST_table.tsv") %>%
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
### ---------------------------------------------------------------------------------------###
### COMBINE
### ---------------------------------------------------------------------------------------###

p_blank <- ggplot() + theme_void()

top_row    <- p_col_ann + p_blank +
  plot_layout(widths = c(10, 0.4))

bottom_row <- p_heat + p_align +
  plot_layout(widths = c(10, 0.4), guides = "collect")

p_final_replicates <- top_row / bottom_row +
  plot_layout(heights = c(0.05, 1)) +
  plot_annotation(
    title    = "Expressed GRC genes — TPM across all replicates",
    subtitle = "Genes clustered by hierarchical clustering | Columns ordered by sex/stage/tissue",
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

p_final_replicates
setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions")
home <- getwd()
ggsave("FigS1_GRC_all_replicate_heatmap_unedited.svg", p_final_replicates, width = 14, height = max(8, nrow(heat_mat) * 0.25 + 3))

