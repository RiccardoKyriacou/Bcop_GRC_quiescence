library(tidyverse)
library(patchwork)

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

# Assign group
df$group <- case_when(
  grepl("TE-like", df$qseqid) ~ "TE",
  grepl("insect-like", df$qseqid) ~ "Insect gene",
  TRUE ~ "Unclear"
)

# Order genes by Alien Index (shared ordering)
df <- df %>%
  arrange(alien_index) %>%
  mutate(gene = factor(gene, levels = gene))

# -----------------------
# LEFT: Alien Index plot
# -----------------------
p1 <- ggplot(
  df,
  aes(
    x = alien_index,
    y = gene,
    colour = interpretation,
    shape = group
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 4) +
  scale_colour_manual(values = c(
    "sciarid-like" = "#41aa96",
    "cecidomyiid-like" = "#aa4196",
    "ambiguous" = "#bdbdbd"
  )) +
  scale_shape_manual(values = c(
    "Insect gene" = 19,
    "TE" = 17,
    "Unclear" = 15
  )) +
  labs(
    x = "Alien Index",
    y = NULL,
    colour = "Ancestry inference",
    shape = "CDS type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "left"
  )

# -----------------------
# RIGHT: Outgroup bitscore barplot
# -----------------------
p2 <- df %>%
  select(gene,
         best_bitscore_cecidomyiid,
         best_bitscore_sciarid,
         best_bitscore_outgroup) %>%
  
  # convert to long format for stacking
  pivot_longer(
    cols = c(best_bitscore_cecidomyiid,
             best_bitscore_sciarid,
             best_bitscore_outgroup),
    names_to = "source",
    values_to = "bitscore"
  ) %>%
  
  mutate(
    source = case_when(
      source == "best_bitscore_cecidomyiid" ~ "Aphidoletes aphidimyza",
      source == "best_bitscore_sciarid" ~ "Bradysia coprophila (core only)",
      source == "best_bitscore_outgroup" ~ "Anopheles gambiae (outgroup)"
    )
  ) %>%
  
  ggplot(aes(
    x = bitscore,
    y = gene,
    fill = source
  )) +
  
  geom_col() +
  
  scale_fill_manual(values = c(
    "Aphidoletes aphidimyza" = "#aa4196",
    "Bradysia coprophila (core only)" = "#41aa96",
    "Anopheles gambiae (outgroup)" = "#333333"
  )) +
  
  labs(
    x = "Best bitscore",
    y = NULL,
    fill = "Reference proteome"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# -----------------------
# COMBINE
# -----------------------
p3 <- p1 + p2 +
  plot_layout(widths = c(2, 1)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

p3

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\Figures\\unedited")
home <- getwd()
home
#ggsave("Alien_Index_plot.svg", p3, width = 10, height = 4.5)

