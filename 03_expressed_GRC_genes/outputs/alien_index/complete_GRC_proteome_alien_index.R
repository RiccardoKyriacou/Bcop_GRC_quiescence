library(tidyverse)
library(patchwork)

# ── Load data ──────────────────────────────────────────────────────────────
df <- read.csv("Full_GRC_proteome_AI_summary.csv", stringsAsFactors = FALSE)
df <- df %>% filter(outgroup_flag != "no_hits_anywhere")

# Shared colour scales
interp_colours <- c(
  "cecidomyiid-like" = "#aa4196",
  "ambiguous"        = "#bdbdbd",
  "sciarid-like"     = "#41aa96"
)

proteome_colours <- c(
  "A. aphidimyza"              = "#aa4196",
  "B. coprophila\n(core only)" = "#41aa96",
  "A. gambiae\n(outgroup)"     = "#333333"
)

df$interpretation <- factor(
  df$interpretation,
  levels = c("cecidomyiid-like", "ambiguous", "sciarid-like")
)

# ── PANEL A: Pie chart ────────────────────────────────────────────────────
pie_counts <- df %>%
  mutate(
    category = case_when(
      interpretation == "cecidomyiid-like" & outgroup_flag == "ok"
      ~ "Cecidomyiid-like",
      interpretation == "cecidomyiid-like" & outgroup_flag == "check_conservation"
      ~ "Cecidomyiid-like\n(conserved across Diptera)",
      interpretation == "ambiguous"
      ~ "Ambiguous",
      interpretation == "sciarid-like" & outgroup_flag == "check_conservation"
      ~ "Sciarid-like\n(conserved across Diptera)",
      interpretation == "sciarid-like" & outgroup_flag == "ok"
      ~ "Sciarid-like"
    )
  ) %>%
  count(category) %>%
  mutate(
    pct   = n / sum(n) * 100,
    label = paste0(formatC(pct, digits = 1, format = "f"), "%\n(n=", n, ")")
  )

pie_levels <- c(
  "Cecidomyiid-like",
  "Cecidomyiid-like\n(conserved across Diptera)",
  "Ambiguous",
  "Sciarid-like\n(conserved across Diptera)",
  "Sciarid-like"
)

pie_pal <- c(
  "Cecidomyiid-like"                             = "#aa4196",
  "Cecidomyiid-like\n(conserved across Diptera)" = "#d48ec8",
  "Ambiguous"                                    = "#bdbdbd",
  "Sciarid-like\n(conserved across Diptera)"     = "#7ecfc5",
  "Sciarid-like"                                 = "#41aa96"
)

pie_counts$category <- factor(pie_counts$category, levels = pie_levels)
pie_counts
pA <- ggplot(pie_counts, aes(x = "", y = n, fill = category)) +
  geom_col(width = 1, colour = "white", linewidth = 0.5) +
  coord_polar(theta = "y", start = 0) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3, colour = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = pie_pal, name = NULL) +
  theme_void(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(1, "cm"),
    legend.direction = "vertical"
  )

pA
# ── PANEL B: AI histogram ─────────────────────────────────────────────────
ai_clip <- 500

pB <- df %>%
  mutate(alien_index = pmax(pmin(alien_index, ai_clip), -ai_clip)) %>%
  ggplot(aes(x = alien_index, fill = interpretation)) +
  geom_histogram(binwidth = 10, colour = NA, alpha = 0.9) +
  geom_vline(xintercept =  0,  linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept =  10, linetype = "dotted", linewidth = 0.4, colour = "grey40") +
  geom_vline(xintercept = -10, linetype = "dotted", linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = interp_colours, name = "Ancestry inference") +
  scale_x_continuous(limits = c(-ai_clip, ai_clip)) +
  annotate("text", x = -380, y = Inf, vjust = 1.8, hjust = 0.5,
           label = "Sciarid-like", colour = "#41aa96",
           size = 3.2, fontface = "italic") +
  annotate("text", x =  380, y = Inf, vjust = 1.8, hjust = 0.5,
           label = "Cecidomyiid-like", colour = "#aa4196",
           size = 3.2, fontface = "italic") +
  labs(x = "Alien Index", y = "No. of GRC transcripts") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
pB

# ── PANEL C: Bitscore violins ─────────────────────────────────────────────
bitscore_long <- df %>%
  select(interpretation,
         best_bitscore_cecidomyiid,
         best_bitscore_sciarid,
         best_bitscore_outgroup) %>%
  pivot_longer(
    cols      = starts_with("best_bitscore"),
    names_to  = "source",
    values_to = "bitscore"
  ) %>%
  mutate(
    source = case_when(
      source == "best_bitscore_cecidomyiid" ~ "A. aphidimyza",
      source == "best_bitscore_sciarid"     ~ "B. coprophila\n(core only)",
      source == "best_bitscore_outgroup"    ~ "A. gambiae\n(outgroup)"
    ),
    source = factor(source, levels = c(
      "A. aphidimyza",
      "B. coprophila\n(core only)",
      "A. gambiae\n(outgroup)"
    ))
  ) %>%
  filter(bitscore > 0) %>%
  mutate(interpretation = factor(
    interpretation,
    levels = c("sciarid-like", "ambiguous", "cecidomyiid-like")
  ))

# Facet labels: capitalise for display
facet_labels <- c(
  "cecidomyiid-like" = "Cecidomyiid-like",
  "ambiguous"        = "Ambiguous",
  "sciarid-like"     = "Sciarid-like"
)

pC <- ggplot(bitscore_long,
             aes(x = source, y = bitscore,
                 fill = source, colour = source)) +
  geom_violin(alpha = 0.6, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               fill = "white", colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ interpretation, nrow = 1, labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = proteome_colours, name = "Reference proteome") +
  scale_colour_manual(values = proteome_colours, guide = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Best bitscore") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(face = "italic", size = 8,
                                    angle = 30, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(face = "italic", size = 8),
    legend.key.size  = unit(0.45, "cm")
  )
pC

# ── COMBINE ───────────────────────────────────────────────────────────────
# A sits left (taller), B and C stack on the right
combined <- pA + (pB / pC) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22))

combined

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\Figures\\unedited")
home <- getwd()
home
ggsave("Full_GRC_AI_combined.svg", combined, width = 13, height = 8)

