library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)
library(dplyr)
library(car)
library(ggtext)

# ---------------------------------------------
# Theme
# ---------------------------------------------

box_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(color = "black", size = 13),
    plot.subtitle = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.title.x  = element_text(color = "black", size = 12),
    axis.text.x   = element_text(size = 11, angle = 45, hjust = 1),
    axis.title.y  = element_text(color = "black", size = 12),
    axis.text.y   = element_text(size = 11),
    legend.position    = "right",
    legend.text        = element_text(size = 11),
    legend.title       = element_text(size = 12),
    strip.text         = element_text(size = 11),
    strip.background   = element_rect(fill = "grey95", color = "grey60"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Colourblind-friendly palette (Wong 2011) — also distinguishable in greyscale
#col_fill  <- c("female" = "#E69F00", "male" = "#56B4E9")
#col_point <- c("female" = "#B37200", "male" = "#1A6E99")

col_fill  <- c("female" = "#F4A5B0", "male" = "#7BB8D4")
col_point <- c("female" = "#C0005C", "male" = "#3B5E8A")

# Data
# ---------------------------------------------

setwd("Revisions\\03_expressed_GRC_genes\\outputs")
home <- getwd()
Expressed_GRC_genes <- read_tsv(file.path(home, "Expressed_GRC_genes.tsv"))

grc_genes <- c("g13362","g13363","g13694","g15244","g17107","g18444",
               "g491","g596","g6396","g7958")

grc_TPM <- Expressed_GRC_genes %>%
  filter(Gene %in% grc_genes)

# ---------------------------------------------
# Shared position object
# ---------------------------------------------

jd <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 42)

# ---------------------------------------------
# Panel A — all expressed GRC genes
# ---------------------------------------------

p0 <- ggplot(Expressed_GRC_genes,
             aes(x = Gene, y = log2_TPM , fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, color = "grey30") +
  geom_point(aes(color = Sex), position = jd, size = 0.5, alpha = 0.7) +
  labs(title = "All expressed GRC genes",
       subtitle = "Only genes with TPM > 0.55 in ≥2/3 of libraries per stage per sex",
       x = "",
       y = expression(TPM~(log[2]))) +
  box_theme +
  theme(legend.position = "top") +
  scale_fill_manual(values = col_fill) +
  scale_color_manual(values = col_point)

# ---------------------------------------------
# Panel B — confidently GRC-linked genes only
# ---------------------------------------------

p1 <- ggplot(grc_TPM,
             aes(x = Gene, y = log2_TPM, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, color = "grey30") +
  geom_point(aes(color = Sex), position = jd, size = 1, alpha = 0.7) +
  labs(title = "Confidently Expressed GRC-linked genes only",
       subtitle = "Only genes with bo somatic mismapping",
       x = "",
       y = expression(TPM~(log[2]))) +
  box_theme +
  theme(legend.position = "top") +
  scale_fill_manual(values = col_fill) +
  scale_color_manual(values = col_point)

# ---------------------------------------------
# Panel C — all genes by developmental stage
# ---------------------------------------------

p2 <- ggplot(Expressed_GRC_genes,
             aes(x = Gene, y = log2_TPM, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, color = "grey30") +
  geom_point(aes(color = Sex), position = jd, size = 0.5, alpha = 0.7) +
  facet_wrap(~Stage, scales = "free_x") +
  labs(title = "All expressed GRC genes — by developmental stage",
       x = "Gene ID",
       y = expression(TPM~(log[2]))) +
  box_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_fill_manual(values = col_fill) +
  scale_color_manual(values = col_point)

# ---------------------------------------------
# Panel D — GRC-linked genes by developmental stage
# ---------------------------------------------

p3 <- ggplot(grc_TPM,
             aes(x = Gene, y = log2_TPM, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, color = "grey30") +
  geom_point(aes(color = Sex), position = jd, size = 1, alpha = 0.7) +
  facet_wrap(~Stage, scales = "free_x") +
  labs(title = "Confidently Expressed GRC genes — by developmental stage",
       x = "Gene ID",
       y = expression(TPM~(log[2]))) +
  box_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_fill_manual(values = col_fill) +
  scale_color_manual(values = col_point)

# ---------------------------------------------
# Combine & save
# ---------------------------------------------

final_plot <- ((p0 | p1) + plot_layout(widths = c(1.5, 1))) /
  ((p2 | p3) + plot_layout(widths = c(1.5, 1))) +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

final_plot

#ggsave("S4_sex_biased.svg", final_plot, width = 13, height = 12)

