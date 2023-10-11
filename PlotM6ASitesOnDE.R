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

limma <- transform(limma, epinano_diff_count=m6a_genes_modified[match(gene_id, m6a_genes_modified$gene_id), ]$epinano_diff_count)
limma[is.na(limma$epinano_diff_count),]$epinano_diff_count <- 0

limma <- transform(limma, xpore_count=m6a_genes_modified[match(gene_id, m6a_genes_modified$gene_id), ]$xpore_count)
limma[is.na(limma$xpore_count),]$xpore_count <- 0

limma_colors <- c("blue", "red", "lightgray")
limma_colors_most_meth <- c("red")

# ------------------ M6ANET

m6a_most_meth_volcano <- function(df, legend_colour, most_methylated, title) {
  ggplot(df, aes(x = logFC, y = neg_log10_P.value, label=gene_id)) +
    geom_point(aes(color = .data[[legend_colour]]), alpha = 1.0) +
    labs(x = "Log Fold Change", y = "-log10 p-value") + 
    theme_classic() +
    scale_colour_manual(values = limma_colors_most_meth) +
    geom_text(aes(label=ifelse(
      match(gene_id, most_methylated) & neg_log10_P.value > 0.5,as.character(gene_id),'')
    ),hjust=0,vjust=0) +
    xlim(-4, 4) +
    ylim(0 , 5) +
    guides(color=guide_legend("")) +
    ggtitle(title)
}

m6a_volcano <- function(df, legend_colour, title) {
  ggplot(df, aes(x = logFC, y = neg_log10_P.value, label=gene_id)) +
    geom_point(aes(color = .data[[legend_colour]]), alpha = 1.0) +
    labs(x = "Log Fold Change", y = "-log10 p-value") + 
    theme_classic() +
    scale_colour_manual(values = limma_colors) +
    xlim(-4, 4) +
    ylim(0 , 5) +
    guides(color=guide_legend("")) +
    ggtitle(title)
}

# ------------------ M6ANET
m6anet_most_methylated_genes <- limma[order(limma$m6anet_count, decreasing=TRUE),]
m6anet_most_methylated_genes <- m6anet_most_methylated_genes[1:10,]
most_methylated_cutoff_m6anet <- min(m6anet_most_methylated_genes$m6anet_count)

limma$num_m6as_in_bin_m6anet <- NA
limma[limma$m6anet_count == 0,]$num_m6as_in_bin_m6anet <- "no methylation"
limma[limma$m6anet_count > 0 & limma$m6anet_count < most_methylated_cutoff_m6anet,]$num_m6as_in_bin_m6anet <- "methylated"
limma[limma$m6anet_count >= most_methylated_cutoff_m6anet,]$num_m6as_in_bin_m6anet <- "most methylated"

m6a_volcano(limma, "num_m6as_in_bin_m6anet", "m6anet")

m6a_most_meth_volcano(
  limma[limma$num_m6as_in_bin_m6anet == "most methylated",], 
  "num_m6as_in_bin_m6anet", 
  m6anet_most_methylated_genes$gene_id,
  "m6anet")

# ------------------ EPINANO
epinano_most_methylated_genes <- limma[order(limma$epinano_count, decreasing=TRUE),]
epinano_most_methylated_genes <- epinano_most_methylated_genes[1:10,]
most_methylated_cutoff_epinano <- min(epinano_most_methylated_genes$epinano_count)

limma$num_m6as_in_bin_epinano <- NA
limma[limma$epinano_count == 0,]$num_m6as_in_bin_epinano <- "no methylation"
limma[limma$epinano_count > 0 & limma$epinano_count < most_methylated_cutoff_epinano,]$num_m6as_in_bin_epinano <- "methylated"
limma[limma$epinano_count >= most_methylated_cutoff_epinano,]$num_m6as_in_bin_epinano <- "most methylated"

m6a_volcano(limma, "num_m6as_in_bin_epinano", "epinano svm")

m6a_most_meth_volcano(
  limma[limma$num_m6as_in_bin_epinano == "most methylated",], 
  "num_m6as_in_bin_epinano", 
  epinano_most_methylated_genes$gene_id,
  "epinano svm")

# ------------------ EPINANO DIFF
epinano_diff_most_methylated_genes <- limma[order(limma$epinano_diff_count, decreasing=TRUE),]
epinano_diff_most_methylated_genes <- epinano_diff_most_methylated_genes[1:10,]
most_methylated_cutoff_epinano_diff <- min(epinano_diff_most_methylated_genes$epinano_diff_count)

limma$num_m6as_in_bin_epinano_diff <- NA
limma[limma$epinano_diff_count == 0,]$num_m6as_in_bin_epinano_diff <- "no methylation"
limma[limma$epinano_diff_count > 0 & limma$epinano_diff_count < most_methylated_cutoff_epinano_diff,]$num_m6as_in_bin_epinano_diff <- "methylated"
limma[limma$epinano_diff_count >= most_methylated_cutoff_epinano_diff,]$num_m6as_in_bin_epinano_diff <- "most methylated"

m6a_volcano(limma, "num_m6as_in_bin_epinano_diff", "epinano diff")

m6a_most_meth_volcano(
  limma[limma$num_m6as_in_bin_epinano_diff == "most methylated",], 
  "num_m6as_in_bin_epinano_diff", 
  epinano_diff_most_methylated_genes$gene_id,
  "epinano diff")

# ------------------ XPORE DIFF
xpore_most_methylated_genes <- limma[order(limma$xpore_count, decreasing=TRUE),]
xpore_most_methylated_genes <- xpore_most_methylated_genes[1:10,]
most_methylated_cutoff_xpore <- min(xpore_most_methylated_genes$xpore_count)

limma$num_m6as_in_bin_xpore <- NA
limma[limma$xpore_count == 0,]$num_m6as_in_bin_xpore <- "no methylation"
limma[limma$xpore_count > 0 & limma$xpore_count < most_methylated_cutoff_xpore,]$num_m6as_in_bin_xpore <- "methylated"
limma[limma$xpore_count >= most_methylated_cutoff_xpore,]$num_m6as_in_bin_xpore <- "most methylated"

m6a_volcano(limma, "num_m6as_in_bin_xpore", "xpore")

m6a_most_meth_volcano(
  limma[limma$num_m6as_in_bin_xpore == "most methylated",], 
  "num_m6as_in_bin_xpore", 
  xpore_most_methylated_genes$gene_id,
  "xpore")
