TPM_genes <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\02_intergenic_TPM\\TPM_genes.tsv")
TPM_genes 

All_GRC_genes <- TPM_genes %>% 
  filter(Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2"))
All_GRC_genes

col_fill  <- c("female" = "#E69F00", "male" = "#56B4E9")
col_point <- c("female" = "#B37200", "male" = "#1A6E99")

MATT_TPM_BH =0.550646
MATT_log2_BH=-0.860803

box_theme <- theme_bw() + 
  theme(
    plot.title = element_text(color = "black", size = 15),
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(color = "black", size = 15),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.title.y = element_text(color = "black", size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    strip.text = element_text(size = 16)
  )



pl0 <- ggplot(All_GRC_genes,
              aes(x = Gene, y = log2_TPM, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, color = "grey30") +
  geom_point(aes(color = Sex),  size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = MATT_log2_BH,
             linetype   = "dashed",
             color      = "red",
             linewidth  = 0.6) +
  annotate("text",
           x     = Inf,
           y     = (MATT_log2_BH + 3),
           label = paste0("TMP cutoff (log2) = ", round(MATT_log2_BH, 2)),
           hjust = 1.05,
           vjust = -0.4,
           size  = 3.2,
           color = "red") +
  labs(title = "TPM values for all GRC-linked genes",
       x = "",
       y = expression(log[2](TPM))) +
  box_theme +
  theme(legend.position = "top") +
  scale_fill_manual(values = col_fill) +
  scale_color_manual(values = col_point)


pl0

