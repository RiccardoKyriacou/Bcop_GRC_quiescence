library(ggplot2)

GRC2_GC <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\long_reads\\GC_content\\GC_SUPER_GRC2_w200000_n100000.tsv")

GRC2_mean_gc=mean(GRC2_GC$GC_Content)

GRC2_GC_mean <- ggplot(GRC2_GC, aes(x = Position, y = GC_Content)) +
  geom_line(color = "skyblue", size = 1.25) +
  geom_hline(yintercept = GRC2_mean_gc, linetype = "dashed", color = "darkblue", size = 1) +
  labs(title = "GRC2",
       subtitle = paste("Mean GC content =", round(GRC2_mean_gc, 2), "%"),
       x = "Position (bp)",
       y = "GC Content (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20, color = "darkblue"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20)
  )
GRC2_GC_mean

flanking_HGT_region_GC <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\long_reads\\GC_content\\GC_651kb_full_HGT_region_w2000_n1000.tsv")

flanking_HGT_mean_gc = mean(flanking_HGT_region_GC$GC_Content)

flanking_HGT_region_GC_mean <- ggplot(flanking_HGT_region_GC, aes(x = Position, y = GC_Content)) +
  geom_line(color = "#4DD0E1", size = 1.25) +
  geom_hline(yintercept = GRC2_mean_gc, linetype = "dashed", color = "darkblue", size = 1) +
  labs(title = "HGT Region + Flanks",
       subtitle = paste("GRC2 GC =", round(GRC2_mean_gc, 2), "%"),
       x = "Position (bp)",
       y = "GC Content (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20, color = "darkblue"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20)
  )
flanking_HGT_region_GC_mean

HGT_region_GC <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\long_reads\\GC_content\\GC_290kb_HGT_region_only_w2000_n1000.tsv")
HGT_mean_gc=mean(HGT_region_GC$GC_Content)

HGT_region_GC_mean <- ggplot(HGT_region_GC, aes(x = Position, y = GC_Content)) +
  geom_line(color = "#00BFC4", size = 1.25) +
  geom_hline(yintercept = HGT_mean_gc, linetype = "dashed", color = "darkblue", size = 1) +
  labs(title = "HGT Region",
       subtitle = paste("Mean GC content =", round(HGT_mean_gc, 2), "%"),
       x = "Position (bp)",
       y = "GC Content (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20, color = "darkblue"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20)
  )
HGT_region_GC_mean


GRC2_GC_mean
GRC2_GC_mean / flanking_HGT_region_GC_mean/ HGT_region_GC_mean
GRC2_GC_mean/HGT_region_GC_mean

### Save 651kb region plot ###
flanking_HGT_region_GC_mean

ggsave("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\figures\\GC_content.svg", 
       plot = flanking_HGT_region_GC_mean, width = 20, height = 7, units = "in", dpi = 300)


### RICKETTSIA GC ###

# Read in the data
rickettsia_GC <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\long_reads\\GC_content\\GC_JAHXDM01.fasta_w2000_n1000.tsv")

scaffold_2 <- rickettsia_GC %>% filter(Scaffold == "ENA|JAHXDM010000002|JAHXDM010000002.1")
scaffold_10 <- rickettsia_GC %>% filter(Scaffold == "ENA|JAHXDM010000010|JAHXDM010000010.1")

mean_2 <- mean(scaffold_2$GC_Content)
mean_10 <- mean(scaffold_10$GC_Content)


plot_2 <- ggplot(scaffold_2, aes(x = Position, y = GC_Content)) +
  geom_line(color = "darkred", size = 1.25) +
  geom_hline(yintercept = mean_2, linetype = "dashed", color = "darkblue", size = 1) +
  labs(
    title = "Scaffold JAHXDM010000002.1",
    subtitle = paste("Mean GC =", round(mean_2, 2), "%"),
    x = "Position (bp)",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16, color = "darkblue"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

plot_10 <- ggplot(scaffold_10, aes(x = Position, y = GC_Content)) +
  geom_line(color = "darkorange3", size = 1.25) +
  geom_hline(yintercept = mean_10, linetype = "dashed", color = "darkblue", size = 1) +
  labs(
    title = "Scaffold JAHXDM010000010.1",
    subtitle = paste("Mean GC =", round(mean_10, 2), "%"),
    x = "Position (bp)",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16, color = "darkblue"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

combined_plot <-HGT_region_GC_mean / plot_2 / plot_10  # vertical layout
combined_plot


