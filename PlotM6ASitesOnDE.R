# Create volcano plot from limma DE analysis
# points should have colour intensity (blue to red) or size (small to big) with 
# respect to how many M6A sites were detected in that gene
library(ggplot2)

limma <- read.csv("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/de/limma-voom_KS-WT.csv", header=TRUE, stringsAsFactors = TRUE)

c1_ks1_plus_epinano_predictions <- read.csv("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/ptmanalysis/epinano_KS1_C1_plus_preds.csv", header=TRUE, stringsAsFactors = TRUE)
c1_ks1_minus_epinano_predictions <- read.csv("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/ptmanalysis/epinano_KS1_C1_minus_preds.csv", header=TRUE, stringsAsFactors = TRUE)
c1_ks1_plus_epinano_predictions <- c1_ks1_plus_epinano_predictions[!duplicated(c1_ks1_plus_epinano_predictions$base_position),]
c1_ks1_minus_epinano_predictions <- c1_ks1_minus_epinano_predictions[!duplicated(c1_ks1_minus_epinano_predictions$base_position),]

c2_m6anet_predictions <- read.csv("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/inference_transcriptome/C2/data.site_proba.csv", header=TRUE, stringsAsFactors = TRUE)
c2_m6anet_predictions <- transform(c2_m6anet_predictions, gene_id=gff[match(transcript_id, gff$ID), ]$gene_id)
m6anet_probabilty_modified_threshold <- 0.8
c2_m6anet_predictions <- c2_m6anet_predictions[c2_m6anet_predictions$probability_modified > m6anet_probabilty_modified_threshold,]


limma$volcano_plot_size <- sample(10, size = nrow(limma), replace = TRUE)
limma$neg_log10_P.value <- -log10(limma$P.Value)

c1_ks1_epinano_predictions <- c1_ks1_plus_epinano_predictions

c1_ks1_epinano_predictions$gene_id <- factor(c1_ks1_epinano_predictions$gene_id, levels=levels(limma$gene_id))

limma$num_m6a <- NA
for(i in 1:nrow(limma)) {
  row <- limma[i,]
  
  matches <- c2_m6anet_predictions[c2_m6anet_predictions$gene_id == row$gene_id,]
  
  limma[i,]$num_m6a <- nrow(matches)
}

limma$num_m6a_norm <- scalar(limma$num_m6a)
limma$num_m6as_in_bin <- NA
limma[limma$num_m6a == 0,]$num_m6as_in_bin <- "0"
limma[limma$num_m6a > 0 & limma$num_m6a <= 2,]$num_m6as_in_bin <- "<= 2"
limma[limma$num_m6a > 2,]$num_m6as_in_bin <- "> 2"

limma_colors <- c("blue", "red", "lightgray") 

ggplot(limma, aes(x = logFC, y = neg_log10_P.value)) +
  geom_point(aes(color = num_m6as_in_bin), alpha = 1.0) +
  labs(x = "Log Fold Change", y = "-log10 p-value") + 
  theme_classic() +
  scale_colour_manual(values = limma_colors)
