library(tidyverse)
library(ggtext)  

# 1. read the ncounts matrices produced by the python script
read_ncounts <- function(file, sp) {
  read_tsv(file, show_col_types = FALSE) %>%
    pivot_longer(-gene, names_to = "chromosome", values_to = "n") %>%
    mutate(species = sp)
}

data <- bind_rows(
  read_ncounts("Bimp_qcov70_pident60_ncounts.tsv", "Bimp"),
  read_ncounts("Ling_qcov70_pident60_ncounts.tsv ", "Ling")
)

# 2. homology class + full gene set (so absent genes still show) 
gene_class <- c(g491="TE", g596="TE", g13363="TE", g17107="TE",
                g6396="gene", g13362="gene", g13694="gene",
                g15244="gene", g18444="gene", g7958="unclear")
all_genes <- names(gene_class)

# order: genes -> unclear -> TEs 
gene_levels <- c("g6396","g13362","g13694","g15244","g18444",
                 "g491","g596","g13363","g17107", 
                 "g7958")

plot_df <- expand_grid(
  gene        = all_genes,
  species     = c("Bimp","Ling"),
  chromosome = c("core_autosome","X","GRC")
) %>%
  left_join(data, by = c("gene","species","chromosome")) %>%
  mutate(n = replace_na(n, 0))

# 3. coloured y-axis labels by class
class_col <- c(TE = "#2c7fb8", gene = "#882255", unclear = "grey40")
y_labels <- setNames(
  sprintf("<span style='color:%s'>%s</span>",
          class_col[gene_class[gene_levels]], gene_levels),
  gene_levels
)

# 4. manual chromosome colours ----
chrm_cols <- c(
  "core_autosome" = "#ffab75",  # orange
  "X"             = "#d65500ff",  # dark orange
  "GRC"           = "#aea0e9ff"   # blue
)
plot_df
# 5. plot 
p3c <- ggplot(filter(plot_df, n > 0),
              aes(chromosome, gene)) +
  geom_point(aes(size = n, fill = chromosome),
             shape = 21, colour = "grey20", stroke = 0.4) +
  facet_wrap(~species,
             labeller = as_labeller(c(Bimp = "B. impatiens",
                                      Ling = "L. ingenua"))) +
  scale_fill_manual(name   = "Chromosome",
                    values = chrm_cols,
                    breaks = c("core_autosome","X","GRC"),
                    labels = c("Autosome","X","GRC")) +
  scale_size_continuous(name = "Loci (n)", range = c(3, 9),
                        breaks = c(1, 2, 4, 8)) +
  scale_x_discrete(limits = c("core_autosome","X","GRC"),
                   labels = c("Core (A)","X","GRC")) +
  scale_y_discrete(limits = rev(gene_levels),
                   labels = y_labels, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_markdown(),
        strip.text  = element_text(face = "italic"))

p3c
ggsave("Fig3C_cross_species.svg", p3c, width = 7, height = 4.5)

