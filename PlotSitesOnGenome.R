# read in gff file to calculate total size of genome bp
# assign resolution variable, ie window size to take avergae of, default should be genome size / 100 for 1% resolution
# count the number of m6a sites in each window size
# plot number of m6a sites against genome position
library(stringr)

# --------------------- CHROMOSOME SUMMARY --------------------- 
write("creating summary of chromosomes... ", stdout())
pfal_gff_path <- "./data/Pfalciparum3D7/gff/data/PlasmoDB-64_Pfalciparum3D7.gff"
rl <- readLines(pfal_gff_path)
gff_chromosome_comments <- rl[grep('^##sequence-region*', rl)]
chromosomes <- data.frame(gff_chromosome_comments)
chromosomes[c("comment", "contig", "start", "end")] <- str_split_fixed(gff_chromosome_comments, " ", 4)
chromosomes$end <- as.numeric(chromosomes$end)
chromosomes$start <- as.numeric(chromosomes$start)

genome_size <- sum(chromosomes$end)

# add genomic start position to chromosomes dfs
chromosomes$genomic_start <- NA
curr_genomic_pos <- 1
for(i in 1:nrow(chromosomes)) {
  chromosomes[i,]$genomic_start <- curr_genomic_pos
  curr_genomic_pos <- curr_genomic_pos + chromosomes[i,]$end
}

chromosomes$genomic_start_percentage <- round(chromosomes$genomic_start / genome_size * 100, 2)

# --------------------- GFF LOAD --------------------- 
# types that have "ID,description"
meta_feature_types <- list("protein_coding_gene", "pseudogene", "ncRNA_gene")
# types that have "ID,parent,..."
feature_types <- list("mRNA", "exon", "CDS", "pseudogenic_transcript", "five_prime_UTR", "three_prime_UTR", "ncRNA", "rRNA", "snoRNA", "tRNA", "snRNA")

# takes in GFF file, splits attributes into gene_ids and description columns
splitAttributes <- function(df) {
  # add new columns: ID, parent, gene_id, attributes
  df[c("ID", "parent", "gene_id", "description")] <- NA
  
  # > unique(pfal_gff$type)
  # [1] "protein_coding_gene"    "mRNA"                   "exon"                  
  # [4] "CDS"                    "pseudogene"             "pseudogenic_transcript"
  # [7] "five_prime_UTR"         "three_prime_UTR"        "ncRNA_gene"            
  # [10] "ncRNA"                  "rRNA"                   "snoRNA"                
  # [13] "tRNA"                   "snRNA"                 
  
  for (feature in meta_feature_types) {
    rows_of_type <- df[df$type == feature,]
    rows_of_type[c("ID", "description")] <- str_split_fixed(rows_of_type$attributes, ";", 2)
    df[df$type == feature,]$ID <- sapply(strsplit(rows_of_type$ID, split='=', fixed=TRUE), function(x) (x[2]))
    df[df$type == feature,]$description <- rows_of_type$description
  }
  
  for (feature in feature_types) {
    rows_of_type <- df[df$type == feature,]
    rows_of_type[c("ID", "parent", "description")] <- str_split_fixed(rows_of_type$attributes, ";", 3)
    df[df$type == feature,]$ID <- sapply(strsplit(rows_of_type$ID, split='=', fixed=TRUE), function(x) (x[2]))
    df[df$type == feature,]$parent <- sapply(strsplit(rows_of_type$parent, split='=', fixed=TRUE), function(x) (x[2]))
    df[df$type == feature,]$description <- rows_of_type$description
  }
  
  # add gene_id to all rows
  # if feature is not parent, gene_id is ID
  df[is.na(df$parent),]$gene_id <- df[is.na(df$parent),]$ID
  # if feature has parent, gene_id is parent without decimal point
  df[!is.na(df$parent),]$gene_id <- sapply(strsplit(df[!is.na(df$parent),]$parent, split='\\.'), function(x) (x[1]))
  
  return(df)
}

write("loading features from gff... ", stdout())
gff <- read.delim(pfal_gff_path, header=FALSE, comment.char="#")
colnames(gff) <- c("contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff <- splitAttributes(gff)
gff <- transform(gff, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
gff$genomic_start_position <- gff$genomic_position + gff$start - 1
gff$length <- gff$end - gff$start
gff$genomic_end_position <- gff$genomic_start_position + gff$length

# --------------------- ADD GENOMIC POSITION --------------------- 

addGenomicPosition <- function(d, type) {
  d <- d
  d$genomic_position <- NA
  
  if (type == "m6anet") {
    # Keep only modified bases above threshold
    m6anet_probabilty_modified_threshold <- 0.8
    d <- d[d$probability_modified > m6anet_probabilty_modified_threshold,]
    
    # clean up data
    d$transcript_start_position <- NA
    d$transcript_end_position <- NA
    d$contig <- NA
    d$strand <- NA
    d <- transform(d, transcript_start_position=gff[match(transcript_id, gff$ID), ]$start)
    d <- transform(d, transcript_end_position=gff[match(transcript_id, gff$ID), ]$end)
    d <- transform(d, contig=gff[match(transcript_id, gff$ID), ]$contig)
    d <- transform(d, strand=gff[match(transcript_id, gff$ID), ]$strand)
    d <- transform(d, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
  
    # NOTE column transcript_position is the start position of the DRACH motiff. it is indexed +1, so transcript_position=37398 maps to contig position=37397
    # with this taken into account, the position of the A in the DRACH motif is only +1 above the column transcript_position
    d$contig_position <- NA
    d[d$strand == "+",]$contig_position <- d[d$strand == "+",]$transcript_start_position+d[d$strand == "+",]$transcript_position
    d[d$strand == "-",]$contig_position <- d[d$strand == "-",]$transcript_end_position-d[d$strand == "-",]$transcript_position
  }
  else if(type == "plasmodb") {
    d <- transform(d, genomic_position=chromosomes[match(Genomic.Sequence.ID, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "epinanosvm") {
    # clean up data
    d$gene_id <- NA
    d[c("start", "end")] <- str_split_fixed(d$Window, "-", 2)
    d$start <- as.numeric(d$start)
    d$end <- as.numeric(d$end)
    # always PLUS 2 from start, epinano doesnt show windows in reverse for - strand like plasmodb
    d$contig_position <- d$start+2
    
    # keep only modified bases
    d <- d[d$prediction == "mod",]
    d <- transform(d, genomic_position=chromosomes[match(Ref, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "epinanodiff") {
    colnames(d) <- c("chr_pos", "ko_feature", "wt_feature", "delta_sum_err", "z_scores", "z_score_prediction")
    d[c("contig", "contig_position", "base", "Strand")] <- str_split_fixed(d$chr_pos, " ", 4)
    d$contig_position <- as.numeric(d$contig_position)
    
    d <- transform(d, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
  }
  else if (type == "xpore") {
    # clean up data
    colnames(d) <- c("id","position","kmer","diff_mod_rate_KO_vs_WT","pval_KO_vs_WT","z_score_KO_vs_WT","mod_rate_KO-rep1","mod_rate_KO-rep2","mod_rate_WT-rep1","mod_rate_WT-rep2","coverage_KO-rep1","coverage_KO-rep2","coverage_WT-rep1","coverage_WT-rep2","mu_unmod","mu_mod","sigma2_unmod","sigma2_mod","conf_mu_unmod","conf_mu_mod","mod_assignment")
  
    d$transcript_start_position <- NA
    d$transcript_end_position <- NA
    d$contig <- NA
    d$strand <- NA
    d <- transform(d, transcript_start_position=gff[match(id, gff$ID), ]$start)
    d <- transform(d, transcript_end_position=gff[match(id, gff$ID), ]$end)
    d <- transform(d, contig=gff[match(id, gff$ID), ]$contig)
    d <- transform(d, strand=gff[match(id, gff$ID), ]$strand)
    d <- transform(d, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
    
    d$contig_position <- NA
    d[d$strand == "+",]$contig_position <- d[d$strand == "+",]$transcript_start_position+d[d$strand == "+",]$position
    d[d$strand == "-",]$contig_position <- d[d$strand == "-",]$transcript_end_position-d[d$strand == "-",]$position
  }

  # add contig position to genomic position and subtract 1.
  # For example in contig 1, position 1: genomic position = 1 + 1 - 1 = 1
  # In contig 2, position 1: genomic position = 640852 + 1 - 1 = 640852
  d$genomic_position <- d$genomic_position + d$contig_position - 1

  return(d)
}

# --------------------- CALCULATE COVERAGE --------------------- 
# count how many DRACH sites there are every sliding_window_size bp genomic_position
num_data_points <- 100
sliding_window_size <- genome_size / num_data_points

calculateCoverage <- function(d, start, end) {
  cur_sliding_window_start <- start
  
  frac_num_points_to_perc <- num_data_points / 100 # 100 percent
  percent_through_genome <- (1:num_data_points) / frac_num_points_to_perc
  gc <- data.frame(percent_through_genome)
  gc$num_sites <- NA
  gc$window_start <- NA
  
  sliding_window_size <- (end-start+1) / num_data_points
  
  for(i in 1:num_data_points) {
    hits <- d[d$genomic_position >= cur_sliding_window_start &
                   d$genomic_position <= cur_sliding_window_start + sliding_window_size,]
    cur_sliding_window_start <- cur_sliding_window_start + sliding_window_size
    
    gc[i,]$num_sites <- nrow(hits)
    gc[i,]$window_start <- cur_sliding_window_start
  }
  
  return(gc)
}

# ---------------------- data from DNA motif search using [GAT][GA]AC[CAT] as pattern in Plasmodium falciparum 3D7
write("calculating DRACH coverage... ", stdout())
drach_positions_path <- "./data/plasmodb/DRACH-DynSpansByMotifSearch_Summary.txt"
plasmodb_drachs <- read.table(file = drach_positions_path, sep = '\t', header = TRUE)
plasmodb_drachs$contig_position <- NA
plasmodb_drachs[plasmodb_drachs$Strand == "-",]$contig_position <- plasmodb_drachs[plasmodb_drachs$Strand == "-",]$Start-2
plasmodb_drachs[plasmodb_drachs$Strand == "+",]$contig_position <- plasmodb_drachs[plasmodb_drachs$Strand == "+",]$Start+2

# add genomic position to plasmodb_drachs
plasmodb_drachs <- addGenomicPosition(plasmodb_drachs, "plasmodb")
drach_genome_site_coverage <- calculateCoverage(plasmodb_drachs, 1, genome_size)

# ---------------------- data from DNA motif search using [GA][GA]AC[CAT] as pattern in Plasmodium falciparum 3D7
write("calculating RRACH coverage... ", stdout())
rrach_positions_path <- "./data/plasmodb/RRACH-DynSpansByMotifSearch_Summary.txt"
plasmodb_rrachs <- read.table(file = rrach_positions_path, sep = '\t', header = TRUE)
plasmodb_rrachs$contig_position <- NA
plasmodb_rrachs[plasmodb_rrachs$Strand == "-",]$contig_position <- plasmodb_rrachs[plasmodb_rrachs$Strand == "-",]$Start-2
plasmodb_rrachs[plasmodb_rrachs$Strand == "+",]$contig_position <- plasmodb_rrachs[plasmodb_rrachs$Strand == "+",]$Start+2

# add genomic position to plasmodb_rrachs
plasmodb_rrachs <- addGenomicPosition(plasmodb_rrachs, "plasmodb")
rrach_genome_site_coverage <- calculateCoverage(plasmodb_rrachs, 1, genome_size)

# ---------------------- M6ANET ----------------------
write("calculating m6anet coverage... ", stdout())

c1_m6anet_path <- "./data/m6anet/C1.data.site_proba.csv"
c2_m6anet_path <- "./data/m6anet/C2.data.site_proba.csv"
ks1_m6anet_path <- "./data/m6anet/KS1.data.site_proba.csv"
ks2_m6anet_path <- "./data/m6anet/KS2.data.site_proba.csv"

c1_m6anet_data <- read.csv(c1_m6anet_path, header=TRUE, stringsAsFactors = TRUE)
c2_m6anet_data <- read.csv(c2_m6anet_path, header=TRUE, stringsAsFactors = TRUE)
ks1_m6anet_data <- read.csv(ks1_m6anet_path, header=TRUE, stringsAsFactors = TRUE)
ks2_m6anet_data <- read.csv(ks2_m6anet_path, header=TRUE, stringsAsFactors = TRUE)

c1_m6anet_data <- addGenomicPosition(c1_m6anet_data, "m6anet")
c2_m6anet_data <- addGenomicPosition(c2_m6anet_data, "m6anet")
ks1_m6anet_data <- addGenomicPosition(ks1_m6anet_data, "m6anet")
ks2_m6anet_data <- addGenomicPosition(ks2_m6anet_data, "m6anet")

c1_m6anet_coverage <- calculateCoverage(c1_m6anet_data, 1, genome_size)
c2_m6anet_coverage <- calculateCoverage(c2_m6anet_data, 1, genome_size)
ks1_m6anet_coverage <- calculateCoverage(ks1_m6anet_data, 1, genome_size)
ks2_m6anet_coverage <- calculateCoverage(ks2_m6anet_data, 1, genome_size)

control_m6anet_coverage <- c1_m6anet_coverage
control_m6anet_coverage$c2_num_sites <- c2_m6anet_coverage$num_sites
control_m6anet_coverage$average_num_sites <- rowMeans(control_m6anet_coverage[,c('num_sites', 'c2_num_sites')], na.rm=TRUE)

ks_m6anet_coverage <- ks1_m6anet_coverage
ks_m6anet_coverage$ks2_num_sites <- ks1_m6anet_coverage$num_sites
ks_m6anet_coverage$average_num_sites <- rowMeans(ks_m6anet_coverage[,c('num_sites', 'ks2_num_sites')], na.rm=TRUE)

# ---------------------- EPINANO SVM ----------------------
write("calculating EpiNano SVM coverage... ", stdout())

# add gene ids to epinano output, only gene ids, not exons etc
just_genes <- gff[is.na(gff$parent),]  
find_gene_id <- function(r) {
  # get all pfal genes that the m6a site is in
  matches <- NA
  matches <- just_genes[
    just_genes$genomic_start_position < as.numeric(getElement(r, "genomic_position"))
    & just_genes$genomic_end_position > as.numeric(getElement(r, "genomic_position"))
    & just_genes$strand == getElement(r, "Strand"),]
  
  if (nrow(matches) == 0) {
    return(NA)
  } 
  else {
    # handles the case where a m6a site is in multiple genes
    return(paste(matches[is.na(matches$parent),]$gene_id, collapse = ','))
  }
}

c1_epinano_path <- "./data/epinano_svm/C1.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
c2_epinano_path <- "./data/epinano_svm/C2.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
ks1_epinano_path <- "./data/epinano_svm/KS1.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
ks2_epinano_path <- "./data/epinano_svm/KS2.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"

c1_epinano_data <- read.csv(c1_epinano_path, header=TRUE, stringsAsFactors = TRUE)
c2_epinano_data <- read.csv(c2_epinano_path, header=TRUE, stringsAsFactors = TRUE)
ks1_epinano_data <- read.csv(ks1_epinano_path, header=TRUE, stringsAsFactors = TRUE)
ks2_epinano_data <- read.csv(ks2_epinano_path, header=TRUE, stringsAsFactors = TRUE)

c1_epinano_data <- addGenomicPosition(c1_epinano_data, "epinanosvm")
c2_epinano_data <- addGenomicPosition(c2_epinano_data, "epinanosvm")
ks1_epinano_data <- addGenomicPosition(ks1_epinano_data, "epinanosvm")
ks2_epinano_data <- addGenomicPosition(ks2_epinano_data, "epinanosvm")

c1_epinano_data$gene_id <- apply(c1_epinano_data, 1, find_gene_id)
c2_epinano_data$gene_id <- apply(c2_epinano_data, 1, find_gene_id)
ks1_epinano_data$gene_id <- apply(ks1_epinano_data, 1, find_gene_id)
ks2_epinano_data$gene_id <- apply(ks2_epinano_data, 1, find_gene_id)

# all m6a sites including ones that couldn't map to genes
c1_epinano_coverage <- calculateCoverage(c1_epinano_data, 1, genome_size)
c2_epinano_coverage <- calculateCoverage(c2_epinano_data, 1, genome_size)
ks1_epinano_coverage <- calculateCoverage(ks1_epinano_data, 1, genome_size)
ks2_epinano_coverage <- calculateCoverage(ks2_epinano_data, 1, genome_size)

control_epinano_coverage <- c1_epinano_coverage
control_epinano_coverage$c2_num_sites <- c2_epinano_coverage$num_sites
control_epinano_coverage$average_num_sites <- rowMeans(control_epinano_coverage[,c('num_sites', 'c2_num_sites')], na.rm=TRUE)

ks_epinano_coverage <- ks1_epinano_coverage
ks_epinano_coverage$ks2_num_sites <- ks2_epinano_coverage$num_sites
ks_epinano_coverage$average_num_sites <- rowMeans(ks_epinano_coverage[,c('num_sites', 'ks2_num_sites')], na.rm=TRUE)

# ---------------------- EPINANO DIFF ----------------------
write("calculating EpiNano Diff coverage... ", stdout())

c1_epinano_diff_path_minus <- "./data/epinano_diff/filtered.KS1_C1_minus_strand.delta-sum_err.prediction.csv"
c1_epinano_diff_path_plus <- "./data/epinano_diff/filtered.KS1_C1_plus_strand.delta-sum_err.prediction.csv"

c2_epinano_diff_path_minus <- "./data/epinano_diff/filtered.KS2_C2_minus_strand.delta-sum_err.prediction.csv"
c2_epinano_diff_path_plus <- "./data/epinano_diff/filtered.KS2_C2_plus_strand.delta-sum_err.prediction.csv"

c1_epinano_diff_data_minus <- read.csv(c1_epinano_diff_path_minus, header=FALSE, stringsAsFactors = TRUE)
c1_epinano_diff_data_plus <- read.csv(c1_epinano_diff_path_plus, header=FALSE, stringsAsFactors = TRUE)
c1_epinano_diff_data <- rbind(c1_epinano_diff_data_minus, c1_epinano_diff_data_plus)

c2_epinano_diff_data_minus <- read.csv(c2_epinano_diff_path_minus, header=FALSE, stringsAsFactors = TRUE)
c2_epinano_diff_data_plus <- read.csv(c2_epinano_diff_path_plus, header=FALSE, stringsAsFactors = TRUE)
c2_epinano_diff_data <- rbind(c2_epinano_diff_data_minus, c2_epinano_diff_data_plus)

c1_epinano_diff_data <- addGenomicPosition(c1_epinano_diff_data, "epinanodiff")
c2_epinano_diff_data <- addGenomicPosition(c2_epinano_diff_data, "epinanodiff")

c1_epinano_diff_data$gene_id <- apply(c1_epinano_diff_data, 1, find_gene_id)
c2_epinano_diff_data$gene_id <- apply(c2_epinano_diff_data, 1, find_gene_id)

c1_epinano_diff_coverage <- calculateCoverage(c1_epinano_diff_data, 1, genome_size)
c2_epinano_diff_coverage <- calculateCoverage(c2_epinano_diff_data, 1, genome_size)

# ---------------------- xPore ----------------------
write("calculating xPore coverage... ", stdout())

xpore_path <- "./data/xpore/majority_direction_kmer_diffmod.table.filtered_pos_diff_mod_rate.DRACHS"
xpore_data <- read.csv(xpore_path, header=FALSE, stringsAsFactors = TRUE)
xpore_data <- addGenomicPosition(xpore_data, "xpore")
xpore_data_coverage <- calculateCoverage(xpore_data, 1, genome_size)

# ---------------------- GRAPH ---------------------- 
write("creating graphs... ", stdout())

out_dir <- "coverage_plots"
dir.create("coverage_plots", showWarnings = FALSE)

plot_width <- 578*2
plot_height <- 458*2
plot_res <- 72*2

# POTENTIAL SITES
jpeg(file = paste(out_dir, "/", "potential_sites_coverage.jpg", sep=""), width = plot_width, height = plot_height, res = plot_res)
plot(
  type = "l",
  drach_genome_site_coverage$percent_through_genome,
  drach_genome_site_coverage$num_sites,
  xlab="% through genome",
  ylab="number of DRACH/RRACH sites",
  ylim=c(0, 8500),
  col="blue"
)
lines(rrach_genome_site_coverage$percent_through_genome, rrach_genome_site_coverage$num_sites, col="red")
for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}
legend(x = "topright",
       legend = c("DRACH", "RRACH"),
       col = c("blue", "red"),
       lwd = 1)
dev.off()

# m6anet
jpeg(file = paste(out_dir, "/", "m6anet_coverage.jpg", sep=""), width = plot_width, height = plot_height, res = plot_res)
plot(
  type = "l",
  control_m6anet_coverage$percent_through_genome,
  control_m6anet_coverage$average_num_sites,
  xlab="% through genome",
  ylab="average number predicted m6a sites",
  col="blue"
)
lines(ks_m6anet_coverage$percent_through_genome, ks_m6anet_coverage$average_num_sites, col="red")

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}
legend(x = "topright",
       legend = c("WT", "KS"),
       col = c("blue", "red"),
       lwd = 1)
dev.off()

# epinano svm
jpeg(file = paste(out_dir, "/", "epinanosvm_coverage.jpg", sep=""), width = plot_width, height = plot_height, res = plot_res)
plot(
  type = "l",
  control_epinano_coverage$percent_through_genome,
  control_epinano_coverage$average_num_sites,
  xlab="% through genome",
  ylab="average number predicted m6a sites",
  col="blue"
)
lines(ks_epinano_coverage$percent_through_genome, ks_epinano_coverage$average_num_sites, col="red")

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}
legend(x = "topright",
       legend = c("WT", "KS"),
       col = c("blue", "red"),
       lwd = 1)
dev.off()

# --------------- COMPARE MODE PLOTS

# epinano diff
jpeg(file = paste(out_dir, "/", "epinano_diff_coverage.jpg", sep=""), width = plot_width, height = plot_height, res = plot_res)
plot(
  type = "l",
  c2_epinano_diff_coverage$percent_through_genome,
  c2_epinano_diff_coverage$num_sites,
  xlab="% through genome",
  ylab="number predicted m6a sites",
  col="red"
)
lines(c1_epinano_diff_coverage$percent_through_genome, c1_epinano_diff_coverage$num_sites, col="blue")

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}
legend(x = "topright",
       legend = c("Rep 1", "Rep 2"),
       col = c("red", "blue"),
       lwd = 1)
dev.off()


# xpore
jpeg(file = paste(out_dir, "/", "xpore_coverage.jpg", sep=""), width = plot_width, height = plot_height, res = plot_res)
plot(
  type = "l",
  xpore_data_coverage$percent_through_genome,
  xpore_data_coverage$num_sites,
  xlab="% through genome",
  ylab="number predicted m6a sites",
  col="blue",
  ylim=c(0,300)
)

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}
legend(x = "topright",
       legend = c("KS-WT replicate predictions"),
       col = c("blue"),
       lwd = 1)
dev.off()

# --------------- get the most modified genes and their numbers
write("determining most modified genes... ", stdout())

# m6anet
m6anet_highlymodified <- sort(table(c2_m6anet_data$transcript_id), decreasing = TRUE)
m6anet_highlymodified <- as.data.frame(m6anet_highlymodified)
colnames(m6anet_highlymodified) <- c("gene_id", "m6anet_count")
m6anet_highlymodified$gene_id <- str_split_fixed(m6anet_highlymodified$gene_id, "\\.", 2)[,1]

# epinano svm
epinano_highlymodified <- sort(table(c2_epinano_data$gene_id), decreasing = TRUE)
epinano_highlymodified <- as.data.frame(epinano_highlymodified)
colnames(epinano_highlymodified) <- c("gene_id", "epinano_count")

# epinano diff
epinano_diff_highlymodified <- sort(table(c1_epinano_diff_data$gene_id), decreasing = TRUE)
epinano_diff_highlymodified <- as.data.frame(epinano_diff_highlymodified)
colnames(epinano_diff_highlymodified) <- c("gene_id", "epinano_diff_count")

# xpore
xpore_highlymodified <- sort(table(xpore_data$id), decreasing = TRUE)
xpore_highlymodified <- as.data.frame(xpore_highlymodified)
colnames(xpore_highlymodified) <- c("gene_id", "xpore_count")
xpore_highlymodified$gene_id <- str_split_fixed(xpore_highlymodified$gene_id, "\\.", 2)[,1]

genes_highlymodified <- merge(m6anet_highlymodified, epinano_highlymodified, by="gene_id")
genes_highlymodified <- merge(genes_highlymodified, epinano_diff_highlymodified, by="gene_id")
genes_highlymodified <- merge(genes_highlymodified, xpore_highlymodified, by="gene_id")

genes_highlymodified$average_count_across_tools <- rowMeans(genes_highlymodified[] [,-1])

m6a_genes_modified <- "./data/m6a_genes_modified.csv"
write.csv(genes_highlymodified, m6a_genes_modified, row.names = FALSE)

# --------------- Pearsons R coeff between treatments on the same tool, based on genome coverage
plot(control_epinano_coverage$average_num_sites, control_m6anet_coverage$average_num_sites)
cor(control_epinano_coverage$average_num_sites, control_m6anet_coverage$average_num_sites)

cor(xpore_data_coverage$num_sites, c1_epinano_diff_coverage$num_sites)
cor(xpore_data_coverage$num_sites, c2_epinano_diff_coverage$num_sites)

# Pearsons R coeff between tools on number of sites detected in each gene
cor(genes_highlymodified$m6anet_count, genes_highlymodified$epinano_count)
cor(rrach_genome_site_coverage$num_sites, drach_genome_site_coverage$num_sites)

# ---------------- get the gene names of the x most highly modified from average across tools
gene_ids_most_mod <- genes_highlymodified[
  order(genes_highlymodified$average_count_across_tools, decreasing = TRUE),
  ][1:20,]$gene_id

eps <- 100
m6anet_most_modified_gene_ids <- genes_highlymodified[
  order(genes_highlymodified$m6anet_count, decreasing = TRUE),
][1:eps,]$gene_id
epinano_diff_most_modified_gene_ids <- genes_highlymodified[
  order(genes_highlymodified$epinano_diff_count, decreasing = TRUE),
][1:eps,]$gene_id
epinano_svm_most_modified_gene_ids <- genes_highlymodified[
  order(genes_highlymodified$epinano_count, decreasing = TRUE),
][1:eps,]$gene_id
xpore_most_modified_gene_ids <- genes_highlymodified[
  order(genes_highlymodified$xpore_count, decreasing = TRUE),
][1:eps,]$gene_id

intersection <- intersect(m6anet_most_modified_gene_ids, epinano_diff_most_modified_gene_ids)
intersection <- intersect(intersection, epinano_svm_most_modified_gene_ids)
intersection <- intersect(intersection, xpore_most_modified_gene_ids)

intersection_without_m6anet <- intersect(epinano_diff_most_modified_gene_ids, epinano_svm_most_modified_gene_ids)
intersection_without_m6anet <- intersect(intersection_without_m6anet, xpore_most_modified_gene_ids)
