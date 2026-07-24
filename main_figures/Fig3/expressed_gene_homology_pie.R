library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# Data
origin_data <- tribble(
  ~Gene,      ~Category,
  "g17107",   "Transposable Element",
  "g13363",   "Transposable Element",
  "g491",     "Transposable Element",
  "g596",     "Transposable Element",
  "g13694",   "Insect",
  "g15244",   "Insect",
  "g7958",    "Unclear",
  "g18444 ",  "Insect",
  "g6396 ",   "Insect",
  "g13362",   "Insect",
)

# Count categories and compute label positions
category_counts <- origin_data %>%
  count(Category) %>%
  arrange(n) %>%
  mutate(
    Category = factor(Category, levels = Category),
    fraction = n / sum(n),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    label_pos = (ymax + ymin) / 2,
    label = paste0(n, "")
  )
category_counts

pie <- ggplot(category_counts, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = Category)) +
  geom_rect(color = "white") +
  geom_text(
    aes(
      x = 1.3,
      y = label_pos,
      label = label
    ),
    size = 5,
    hjust = 0.5
  ) +
  coord_polar(theta = "y") +
  xlim(c(-0.5, 1.6)) +
  theme_void() +
  scale_fill_manual(
    values = c(
      "Transposable Element" = "#88cceeff",
      "Insect" = "#882255ff",
      "Uncertain" = "grey65"
    )
  ) +
  labs(fill = "Top BLAST Hit") +
  theme(legend.position = "bottom")

pie

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\03_expressed_GRC_genes\\outputs")
home <- getwd()
ggsave("GRC_BLAST_pie.svg", pie,
       width = 14, height = 10)
