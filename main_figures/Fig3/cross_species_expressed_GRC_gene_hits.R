
library(tidyverse)

file <- "outputs\\Bimp_Ling_BLAST\\Bimp_ling_combined_blast_table.tsv"

blast <- read.delim(file, header = FALSE)

colnames(blast) <- c(
  "gene","subject","pident","length","mismatch","gapopen",
  "qstart","qend","sstart","send","evalue","bitscore",
  "qlen","slen"
)
blast
# -----------------------
# P4: Cross-species conservation
# -----------------------
blast

# Classify hits
blast_scaffolds <- blast %>%
  mutate(
    category = case_when(
      grepl("^Bimp_", subject) & grepl("GRC", subject) ~ "B. impatiens GRC",
      grepl("^Bimp_", subject) ~ "B. impatiens Core",
      grepl("^Ling_", subject) & grepl("GRC", subject) ~ "L. ingenua GRC",
      grepl("^Ling_", subject) ~ "L. ingenua Core",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  
  # Count unique chromosomes/scaffolds
  distinct(gene, category, subject) %>%
  count(gene, category, name = "n_scaffolds")


# Add genes with no hits

all_genes <- unique(df$gene)


blast_scaffolds <- blast_scaffolds %>%
  complete(
    gene = all_genes,
    category = c(
      "B. impatiens GRC",
      "B. impatiens Core",
      "L. ingenua GRC",
      "L. ingenua Core"
    ),
    fill = list(n_scaffolds = 0)
  )
blast_scaffolds

# Identify genes with absolutely no hits
genes_with_hits <- blast_scaffolds %>%
  group_by(gene) %>%
  summarise(total_hits = sum(n_scaffolds), .groups = "drop")

no_hit_genes <- genes_with_hits %>%
  filter(total_hits == 0) %>%
  transmute(
    gene,
    category = "No hit",
    n_scaffolds = 1
  )

blast_scaffolds <- bind_rows(
  blast_scaffolds,
  no_hit_genes
)

# Match ordering from p1

# Read table
df <- read.delim(
  "C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\03_expressed_GRC_genes\\outputs\\alien_index\\Alien_Index_summary.txt"
)

# Clean gene names
df$gene <- gsub(
  "_(insect-like|TE-like|unclear_homology)",
  "",
  df$qseqid
)


gene_order <- df %>%
  distinct(gene, alien_index) %>%
  arrange(alien_index) %>%
  pull(gene)

blast_scaffolds$gene <- factor(
  blast_scaffolds$gene,
  levels = gene_order
)

p4 <- ggplot(
  blast_scaffolds,
  aes(
    x = n_scaffolds,
    y = gene,
    fill = category
  )
) +
  geom_col() +
  
  scale_fill_manual(
    values = c(
      "B. impatiens GRC" = "#aea0e9ff",
      "B. impatiens Core" = "#ffab75ff",
      "L. ingenua GRC1" = "#321e87ff",
      "L. ingenua Core" = "#d65500ff",
      "No hit" = "white"
    )
  ) +
  
  scale_y_discrete(drop = FALSE) +
  
  labs(
    x = "No. of cross-species hits",
    y = NULL,
    fill = "Cross-species conservation"
  ) +
  
  theme_classic(base_size = 14)
p4

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\Figures\\unedited")
home <- getwd()
home
ggsave("Cross_species_hits.svg", p4, width = 5, height = 4)

