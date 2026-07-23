# Load required libraries
library(mclust)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)

### ---------------------------------------------------------------------------------------###
### STEP 1) Load and prepare data ###
### ---------------------------------------------------------------------------------------###
#file.choose()
TPM_intergenic <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\02_intergenic_TPM\\outputs\\combined_intergenic_TPM.tsv", 
                           col_names = FALSE)
colnames(TPM_intergenic) <- c("Chromosome", "TPM", "Gene", "Species", "Sex", "Tissue", "Stage", "Sample")
TPM_intergenic$TPM <- as.numeric(TPM_intergenic$TPM)
# Take log and add a small offset to avoid log2(0) issues
TPM_intergenic$log2_TPM <- log2(TPM_intergenic$TPM + 1e-6)
TPM_intergenic

### ---------------------------------------------------------------------------------------###
### STEP 2) Perform Gaussian mixture model (GMM) fitting ###
### ---------------------------------------------------------------------------------------###
set.seed(33)  # Ensures reproducible EM initialisation
model <- Mclust(TPM_intergenic$log2_TPM, G = 3)
#model <- Mclust(TPM_intergenic$log2_TPM)
summary(model)
plot(model, what = "BIC")  # Plot BIC values for different component numbers
# Extract means, standard deviations, and proportions of each component
means <- model$parameters$mean
sds <- sqrt(model$parameters$variance$sigma)  
proportions <- model$parameters$pro  

# Assign each data point to a Gaussian component
TPM_intergenic$Component <- model$classification  

### ---------------------------------------------------------------------------------------###
### STEP 3) Define species-specific maximum inactively transcribed threshold (MITTspecies) ###
### ---------------------------------------------------------------------------------------###

# Find the component with lowest overlap to CDS curve (assumed to be non-expressed regions)
lowest_component_overlap <- 1

# Find the highest TPM value among points assigned to this component
MITTspecies_log2 <- max(TPM_intergenic$log2_TPM[TPM_intergenic$Component == lowest_component_overlap])
MITTspecies_log2
# e.g. 99th percentile of Component 1
#MITTspecies_log2 <- quantile(TPM_intergenic$log2_TPM[TPM_intergenic$Component == 1], 0.99)
#MITTspecies_log2

MITTspecies <- 2^MITTspecies_log2  # Convert back to TPM scale

print(MITTspecies)  # This is the final threshold

### ---------------------------------------------------------------------------------------###
### STEP 4) PLOT  ###
### ---------------------------------------------------------------------------------------###

# Prepare density data
dens <- density(TPM_intergenic$log2_TPM)

# Create a dataframe for Gaussian components
x_vals <- seq(min(dens$x), max(dens$x), length.out = 1000)
gaussians <- data.frame(
  x = rep(x_vals, length(means)),
  component = factor(rep(1:length(means), each = length(x_vals))),
  y = unlist(lapply(1:length(means), function(i) {
    dnorm(x_vals, mean = means[i], sd = sds[i]) * proportions[i]
  }))
)

# Generate the plot
gaussian_plot <- ggplot(data = TPM_intergenic, aes(x = log2_TPM)) +
  # Plot density of intergenic TPM values
  geom_density(color = "black", fill = "salmon", alpha = 0.3) +
  
  # Add Gaussian components as separate lines
  geom_line(data = gaussians, aes(x = x, y = y, color = component), size = 1.2) +
  scale_color_manual(values = rainbow(length(means))) + 
  #scale_color_manual(values = scales::viridis_pal(option = "inferno")(length(means))) +  # different colour potion
  
  # Add vertical threshold line for MITTspecies
  geom_vline(xintercept = log2(MITTspecies), color = "red", linetype = "dashed", size = 1)+
  labs(title = "Intergenic Log2(TPM) Density with Gaussian Components",
       x = "Log2(TPM)",
       y = "Density",
       color = "Gaussian") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right") +
  annotate("text", x = log2(MITTspecies) + 0.5, y = 0.22,
           label = paste("MITTspecies =", round(log2(MITTspecies), 2)),
           color = "red", hjust = 0)

gaussian_plot

### ---------------------------------------------------------------------------------------###
### STEP 5) Filter intergenic regions below MITTspecies threshold
### ---------------------------------------------------------------------------------------###

reference_intergenic <- TPM_intergenic[TPM_intergenic$log2_TPM < log2(MITTspecies), ]
head(reference_intergenic)

# Compute mean and standard deviation of reference intergenic regions
mean_ref <- mean(reference_intergenic$log2_TPM, na.rm = TRUE)
sd_ref <- sd(reference_intergenic$log2_TPM, na.rm = TRUE)

### ---------------------------------------------------------------------------------------###
### STEP 6) Load gene TPM data and compute Z-scores for genes
### ---------------------------------------------------------------------------------------###
# Load gene expression data 
TPM_genes <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\01_STAR_TPM\\combined_TPM_only.tsv", 
                      col_names = FALSE)
colnames(TPM_genes) <- c("Chromosome", "TPM", "Gene", "Species", "Sex", "Tissue", "Stage", "Sample")
TPM_genes$TPM <- as.numeric(TPM_genes$TPM)
#Apply correction to 0 values
TPM_genes$log2_TPM <- log2(TPM_genes$TPM + 1e-6)
TPM_genes
### ---------------------------------------------------------------------------------------###
### Step 7) Compute Z-scores for genes
### ---------------------------------------------------------------------------------------###
# Compute Z-score
TPM_genes$Z_score <- (TPM_genes$log2_TPM - mean_ref) / sd_ref

### ---------------------------------------------------------------------------------------###
### Step 8) Compute -pvalues for genes
### ---------------------------------------------------------------------------------------###
# Compute one-sided p-value
TPM_genes$p_value <- 1 - pnorm(TPM_genes$Z_score)
TPM_genes

### ---------------------------------------------------------------------------------------###
### Step 9) Apply Benjamini-Hochberg (BH) correction
### ---------------------------------------------------------------------------------------###
TPM_genes <- TPM_genes %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

### ---------------------------------------------------------------------------------------###
### Step 10) # Identify MATT_TPM_BH threshold (smallest TPM where BH-corrected p-value ≤ alpha)
### ---------------------------------------------------------------------------------------###
alpha <- 0.05  # Define -pvalue threshold (p = 0.05)

MATT_log2_BH <- min(TPM_genes$log2_TPM[TPM_genes$p_adj <= alpha], na.rm = TRUE)
# Convert back to TPM scale
MATT_TPM_BH <- 2^MATT_log2_BH

# Classify genes using MATT
TPM_genes <- TPM_genes %>%
  mutate(Expression_Status_BH = ifelse(TPM >= MATT_TPM_BH, "Active", "Inactive"))
TPM_genes
# Output final threshold
MATT_TPM_BH
MATT_log2_BH
### ---------------------------------------------------------------------------------------###
### Step 11) Plot gene expression vs. reference intergenic
### ---------------------------------------------------------------------------------------###

# Classify genes as "Expressed" or "Not Expressed" based on MATT_TPM_BH
TPM_genes$Expression_Status <- ifelse(TPM_genes$TPM >= MATT_TPM_BH, "Expressed", "Not Expressed")
table(TPM_genes$Expression_Status)

p_classified <- ggplot(TPM_genes, aes(x = log2_TPM, fill = Expression_Status)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("skyblue", "darkblue")) +  
  geom_vline(xintercept = log2(MATT_TPM_BH), color = "black", linetype = "dashed", size = 1) +
  annotate("text", x = log2(MATT_TPM_BH) + 0.5, y = 25000,
           label = paste("MATTlibrary =", round(log2(MATT_TPM_BH), 2)),
           color = "black", hjust = 0)+
  labs(title = "Gene Expression Classification Based on MATTlibrary",
       x = "Log2(TPM)",
       y = "Count",
       fill = "Expression Status") +
  theme_minimal()

print(p_classified)

# For GRC
TPM_genes_GRC <- TPM_genes %>% filter(Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2"))
TPM_genes_GRC

# Plot only GRC genes
p_classified_GRC <- ggplot(TPM_genes_GRC, aes(x = log2_TPM, fill = Expression_Status)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("#86E7B8", "darkblue")) +  
  annotate("text", x = log2(MATT_TPM_BH) + 0.5, y = 45,
           label = paste("MATTlibrary =", round(log2(MATT_TPM_BH), 2)),
           color = "black", hjust = 0)+
  geom_vline(xintercept = log2(MATT_TPM_BH), color = "black", linetype = "dashed", size = 1) +
  labs(title = "GRC Gene Expression Classification Based on MATTlibrary",
       x = "Log2(TPM)",
       y = "Count",
       fill = "Expression Status") +
  theme_minimal()
print(p_classified_GRC)

# Count genes on each scaffold
table(TPM_genes$Chromosome[TPM_genes$Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2")])
# Count expressed genes on each scaffold
table(TPM_genes$Chromosome[TPM_genes$Chromosome %in% c("SUPER_GRC1", "SUPER_GRC2") & TPM_genes$Expression_Status == "Expressed"])

### ---------------------------------------------------------------------------------------###
### Step 12) Full plot
### ---------------------------------------------------------------------------------------###

#Add type column to data
TPM_genes$Type <- "Gene"  
mean_gene_TPM <- mean(TPM_genes$log2_TPM, na.rm = TRUE)
TPM_intergenic$Type <- "Intergenic"
mean_intergenic_TPM <- mean(TPM_intergenic$log2_TPM, na.rm = TRUE)

# Combine both datasets
combined_TPM_data <- bind_rows(TPM_genes, TPM_intergenic)
combined_TPM_data

# Density plot with both distributions on one graph
p_combined <- ggplot(combined_TPM_data, aes(x = log2(TPM), fill = Type)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = mean_gene_TPM, color = "darkblue", linetype = "dashed", size = 1) + 
  annotate("text", x = mean_gene_TPM, y = 0.2,
           label = paste("x̄ =", round(mean_gene_TPM, 2)),
           color = "darkblue", hjust = -0.1) +
  geom_vline(xintercept = mean_intergenic_TPM, color = "red", linetype = "dashed", size = 1) + 
  annotate("text", x = mean_intergenic_TPM, y = 0.2,
           label = paste("x̄ =", round(mean_intergenic_TPM, 2)),
           color = "red", hjust = -0.1) +
  labs(title = "Genic vs Intergenic TPM distributions", x = "Log2(TPM)", y = "Density") +
  scale_fill_manual(values = c("skyblue", "salmon")) +  
  theme_minimal()
p_combined

### Plot Theme ###
box_theme <- theme_bw() + 
  theme(plot.title = element_text(color="black", size=15),
        panel.background = element_rect(fill="white"),
        axis.title.x = element_text(color="black", size=15),
        axis.text.x = element_text(size = 15, hjust = 1),
        axis.title.y = element_text(color="black", size=15),
        axis.text.y = element_text(size = 15), 
        legend.position = "right", 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16)
  )

# Apply the same x and y limits to both plots
gaussian_plot <- gaussian_plot + 
  xlim(-5, 15) + box_theme
p_combined <- p_combined + 
  xlim(-5, 15) + box_theme 
p_classified <- p_classified + 
  xlim(-5, 15) + box_theme
p_classified_GRC <- p_classified_GRC + 
  xlim(-5, 15) + box_theme

# Combine both plots using patchwork ###
plot_final_bcop <- p_combined + gaussian_plot + p_classified + p_classified_GRC
plot_final_bcop

### ---------------------------------------------------------------------------------------###
### Step 13) Full plot with CDS outlines
### ---------------------------------------------------------------------------------------###

# Prepare component names to match legend labels
gaussians$comp_name <- paste0("Component ", gaussians$component)

# Create named color vector: one for CDS and one per component
comp_names <- unique(gaussians$comp_name)
color_values <- c("CDS Density" = "deepskyblue2",
                  setNames(rainbow(length(comp_names)), comp_names))

# Create named linetype vector: dashed for CDS, solid for components
linetype_values <- c("CDS Density" = "dashed",
                     setNames(rep("solid", length(comp_names)), comp_names))

# Compute density for gene log2(TPM) values
gene_dens <- density(TPM_genes$log2_TPM)

# Prepare a data frame for the gene density curve
gene_density_df <- data.frame(x = gene_dens$x, y = gene_dens$y, type = "CDS Density")

# Plot gene density and intergenic Gaussian components together (with legend entry for CDS)
P_cds <- ggplot() +
  # Gene density curve mapped to a name so it appears in the legend
  geom_line(data = gene_density_df,
            aes(x = x, y = y, color = "CDS Density", linetype = "CDS Density"),
            size = 1.3) +
  # Intergenic density curve (filled salmon) — no legend entry (remove mapping)
  geom_density(data = TPM_intergenic, aes(x = log2_TPM),
               fill = "salmon", alpha = 0.3, color = "black", size = 0.8) +
  # Add each Gaussian component as a separate line, mapped to "Component N" names
  geom_line(data = gaussians, aes(x = x, y = y, color = comp_name, linetype = comp_name), size = 1.2) +
  # Manual scales — ensure names align with the mapped labels above
  scale_color_manual(name = "Curve Type", values = color_values) +
  scale_linetype_manual(name = "Curve Type", values = linetype_values) +
  labs(title = "Overlay: Gene Density & Intergenic Gaussian Components",
       x = "Log2(TPM)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "right") + box_theme +
  annotate("text", x = log2(MITTspecies) + 0.5, y = 0.22, 
           label = paste("MITTspecies =", round(log2(MITTspecies), 2)),
           color = "red", hjust = 0) + 
  geom_vline(xintercept = log2(MITTspecies), color = "red", linetype = "dashed", size = 1)

P_cds

################
# FINAL FIGURE
################

plot_final_bcop <- p_combined + P_cds + p_classified + p_classified_GRC

plot_final_bcop

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\02_intergenic_TPM")
home <- getwd()
home
#ggsave(plot = plot_final_bcop, filename = file.path(home, "S1_Background_unedited_TPM.svg"), 
       #device = "svg", units = "cm", width = 45, height = 30)

 
################
# Save output
################)
TPM_genes
write.table(TPM_genes,
            file = file.path(home, "TPM_genes.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

TPM_genes



