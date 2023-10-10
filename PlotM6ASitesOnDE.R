# Create volcano plot from limma DE analysis
# points should have colour intensity (blue to red) or size (small to big) with 
# respect to how many M6A sites were detected in that gene
library(ggplot2)

limma <- read.csv("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/de/limma-voom_KS-WT.csv", header=TRUE, stringsAsFactors = TRUE)
m6a_genes_modified_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6a_genes_modified.csv"

m6a_genes_modified <- read.csv(m6a_genes_modified_path, header=TRUE, stringsAsFactors = TRUE)

limma$volcano_plot_size <- sample(10, size = nrow(limma), replace = TRUE)
limma$neg_log10_P.value <- -log10(limma$P.Value)

m6a_genes_modified$gene_id <- factor(m6a_genes_modified$gene_id, levels=levels(limma$gene_id))

limma <- transform(limma, m6anet_count=m6a_genes_modified[match(gene_id, m6a_genes_modified$gene_id), ]$m6anet_count)
limma[is.na(limma$m6anet_count),]$m6anet_count <- 0

limma <- transform(limma, epinano_count=m6a_genes_modified[match(gene_id, m6a_genes_modified$gene_id), ]$epinano_count)
limma[is.na(limma$epinano_count),]$epinano_count <- 0

# ------------------ M6ANET
m6anet_most_methylated_genes <- limma[order(limma$m6anet_count, decreasing=TRUE),]
m6anet_most_methylated_genes <- m6anet_most_methylated_genes[1:10,]
most_methylated_cutoff_m6anet <- min(m6anet_most_methylated_genes$m6anet_count)

limma$num_m6as_in_bin_m6anet <- NA
limma[limma$m6anet_count == 0,]$num_m6as_in_bin_m6anet <- "no methylation"
limma[limma$m6anet_count > 0 & limma$m6anet_count <= most_methylated_cutoff_m6anet,]$num_m6as_in_bin_m6anet <- "methylated"
limma[limma$m6anet_count > most_methylated_cutoff_m6anet,]$num_m6as_in_bin_m6anet <- "most methylated"


limma_colors <- c("blue", "red", "lightgray")

ggplot(limma, aes(x = logFC, y = neg_log10_P.value, label=gene_id)) +
  geom_point(aes(color = num_m6as_in_bin_m6anet), alpha = 0.8) +
  labs(x = "Log Fold Change", y = "-log10 p-value") + 
  theme_classic() +
  scale_colour_manual(values = limma_colors) +
  geom_text(aes(label=ifelse(match(gene_id, m6anet_most_methylated_genes$gene_id),as.character(gene_id),'')),hjust=0,vjust=0)

# ------------------ EPINANO
epinano_most_methylated_genes <- limma[order(limma$epinano_count, decreasing=TRUE),]
epinano_most_methylated_genes <- epinano_most_methylated_genes[1:10,]
most_methylated_cutoff_epinano <- min(epinano_most_methylated_genes$epinano_count)

limma$num_m6as_in_bin_epinano <- NA
limma[limma$epinano_count == 0,]$num_m6as_in_bin_epinano <- "no methylation"
limma[limma$epinano_count > 0 & limma$epinano_count <= most_methylated_cutoff_epinano,]$num_m6as_in_bin_epinano <- "methylated"
limma[limma$epinano_count > most_methylated_cutoff_epinano,]$num_m6as_in_bin_epinano <- "most methylated"

limma_colors <- c("blue", "red", "lightgray")

ggplot(limma, aes(x = logFC, y = neg_log10_P.value, label=gene_id)) +
  geom_point(aes(color = num_m6as_in_bin_epinano), alpha = 0.8) +
  labs(x = "Log Fold Change", y = "-log10 p-value") + 
  theme_classic() +
  scale_colour_manual(values = limma_colors) +
  geom_text(aes(label=ifelse(match(gene_id, epinano_most_methylated_genes$gene_id),as.character(gene_id),'')),hjust=0,vjust=0)
