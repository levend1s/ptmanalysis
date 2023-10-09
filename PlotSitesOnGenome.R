# read in gff file to calculate total size of genome bp
# assign resolution variable, ie window size to take avergae of, default should be genome size / 100 for 1% resolution
# count the number of m6a sites in each window size
# plot number of m6a sites against genome position

pfal_gff_path <- "~/Desktop/Biomedical\ Research\ Project/Pfalciparum3D7/gff/data/PlasmoDB-64_Pfalciparum3D7.gff"
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

addGenomicPosition <- function(d, type) {
  d <- d
  d$genomic_position <- NA
  
  if (type == "m6anet") {
    d <- transform(d, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "plasmodb") {
    d <- transform(d, genomic_position=chromosomes[match(Genomic.Sequence.ID, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "epinanosvm") {
    d <- transform(d, genomic_position=chromosomes[match(Ref, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "epinano") {
    d <- transform(d, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
  }

  # add contig position to genomic position and subtract 1. This handles the case where genomic start = 1 and pos = 1 -> 1+1 - 1.
  d$genomic_position <- d$genomic_position + d$pos - 1

  return(d)
}

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

gff <- read.delim(pfal_gff_path, header=FALSE, comment.char="#")
colnames(gff) <- c("contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# types that have "ID,description"
meta_feature_types <- list("protein_coding_gene", "pseudogene", "ncRNA_gene")
# types that have "ID,parent,..."
feature_types <- list("mRNA", "exon", "CDS", "pseudogenic_transcript", "five_prime_UTR", "three_prime_UTR", "ncRNA", "rRNA", "snoRNA", "tRNA", "snRNA")

gff <- splitAttributes(gff)

# count how many DRACH sites there are every sliding_window_size bp genomic_position
num_data_points <- 100
sliding_window_size <- genome_size / num_data_points
frac_num_points_to_perc <- num_data_points / 100
percent_through_genome <- (1:num_data_points) / frac_num_points_to_perc

calculateCoverage <- function(d) {
  cur_sliding_window_start <- 1
  gc <- data.frame(percent_through_genome)
  gc$num_sites <- NA
  gc$window_start <- NA
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
drach_positions_path <- "/Users/joshualevendis/Downloads/DRACH-DynSpansByMotifSearch_Summary.txt"
plasmodb_drachs <- read.table(file = drach_positions_path, sep = '\t', header = TRUE)
plasmodb_drachs$pos <- NA
plasmodb_drachs[plasmodb_drachs$Strand == "-",]$pos <- plasmodb_drachs[plasmodb_drachs$Strand == "-",]$Start-2
plasmodb_drachs[plasmodb_drachs$Strand == "+",]$pos <- plasmodb_drachs[plasmodb_drachs$Strand == "+",]$Start+2

# add genomic position to plasmodb_drachs
plasmodb_drachs <- addGenomicPosition(plasmodb_drachs, "plasmodb")

drach_genome_site_coverage <- calculateCoverage(plasmodb_drachs)

plot(
  type = "l",
  drach_genome_site_coverage$percent_through_genome,
  drach_genome_site_coverage$num_sites,
  xlab="% through genome",
  ylab="number of DRACH sites"
)
for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}

# ---------------------- data from DNA motif search using [GA][GA]AC[CAT] as pattern in Plasmodium falciparum 3D7
rrach_positions_path <- "/Users/joshualevendis/Downloads/RRACH-DynSpansByMotifSearch_Summary.txt"
plasmodb_rrachs <- read.table(file = rrach_positions_path, sep = '\t', header = TRUE)
plasmodb_rrachs$pos <- NA
plasmodb_rrachs[plasmodb_rrachs$Strand == "-",]$pos <- plasmodb_rrachs[plasmodb_rrachs$Strand == "-",]$Start-2
plasmodb_rrachs[plasmodb_rrachs$Strand == "+",]$pos <- plasmodb_rrachs[plasmodb_rrachs$Strand == "+",]$Start+2

# add genomic position to plasmodb_rrachs
plasmodb_rrachs <- addGenomicPosition(plasmodb_rrachs, "plasmodb")

rrach_genome_site_coverage <- calculateCoverage(plasmodb_rrachs)

plot(
  type = "l",
  rrach_genome_site_coverage$percent_through_genome,
  rrach_genome_site_coverage$num_sites,
  xlab="% through genome",
  ylab="number of RRACH sites"
)
for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}

# ---------------------- M6ANET ----------------------
generateM6AnetCoverage <- function(path) {
  #path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/inference_transcriptome/C2/data.site_proba.csv"
  m6anet_data <- read.csv(path, header=TRUE, stringsAsFactors = TRUE)
  # NOTE column transcript_position is the start position of the DRACH motiff. it is indexed +1, so transcript_position=37398 maps to contig position=37397
  # with this taken into account, the position of the A in the DRACH motif is only +1 above the column transcript_position
  m6anet_data$transcript_start_position <- NA
  m6anet_data$contig <- NA
  m6anet_data <- transform(m6anet_data, transcript_start_position=gff[match(transcript_id, gff$ID), ]$start)
  m6anet_data$position <- m6anet_data$transcript_position+m6anet_data$transcript_start_position+1
  m6anet_data <- transform(m6anet_data, contig=gff[match(transcript_id, gff$ID), ]$contig)
  
  m6anet_probabilty_modified_threshold <- 0.8
  m6anet_data_filtered <- m6anet_data[m6anet_data$probability_modified > m6anet_probabilty_modified_threshold,]
  
  # add genomic position
  m6anet_data_filtered <- addGenomicPosition(m6anet_data_filtered, "m6anet")
  m6anet_genome_site_coverage <- calculateCoverage(m6anet_data_filtered)
  
  return(m6anet_genome_site_coverage)
}

#c1_m6anet_coverage <- generateM6AnetCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/c1_m6anet_output_inference/c1_data.site_proba.csv")
c2_m6anet_coverage <- generateM6AnetCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/inference_transcriptome/C2/data.site_proba.csv")
#ks1_m6anet_coverage <- generateM6AnetCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/ks_1_m6anet_output_inference/ks1_data.site_proba.csv")
ks2_m6anet_coverage <- generateM6AnetCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/inference_transcriptome/KS2/data.site_proba.csv")

control_m6anet_coverage <- c2_m6anet_coverage
#control_m6anet_coverage$c2_num_sites <- c2_m6anet_coverage$num_sites
#control_m6anet_coverage$average_num_sites <- rowMeans(control_m6anet_coverage[,c('num_sites', 'c2_num_sites')], na.rm=TRUE)

ks_m6anet_coverage <- ks2_m6anet_coverage
#ks_m6anet_coverage$ks2_num_sites <- ks2_m6anet_coverage$num_sites
#ks_m6anet_coverage$average_num_sites <- rowMeans(ks_m6anet_coverage[,c('num_sites', 'ks2_num_sites')], na.rm=TRUE)

# GRAPH
plot(
  type = "l",
  control_m6anet_coverage$percent_through_genome,
  control_m6anet_coverage$num_sites,
  xlab="% through genome",
  ylab="number of DRACH sites",
  col="blue"
)
lines(ks_m6anet_coverage$percent_through_genome, ks_m6anet_coverage$num_sites, col="red")

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}

# ---------------------- EPINANO ----------------------
generateEpinanoCoverage <- function(path) {
  epinano_data <- read.csv(path, header=TRUE, stringsAsFactors = TRUE)
  epinano_data[c("start", "end")] <- str_split_fixed(epinano_data$Window, "-", 2)
  epinano_data$start <- as.numeric(epinano_data$start)
  epinano_data$end <- as.numeric(epinano_data$end)
  # always PLUS 2 from start, epinano doesnt show windows in reverse for - strand like plasmodb
  epinano_data$position <- epinano_data$start+2
  
  epinano_data <- epinano_data[epinano_data$prediction == "mod",]
  
  epinano_data <- addGenomicPosition(epinano_data, "epinanosvm")
  epinano_genome_site_coverage <- calculateCoverage(epinano_data)
  
  return(epinano_genome_site_coverage)
}

c1_epinano_coverage <- generateEpinanoCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/C1.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv")
c2_epinano_coverage <- generateEpinanoCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/C2.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv")
ks1_epinano_coverage <- generateEpinanoCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/KS1.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv")
ks2_epinano_coverage <- generateEpinanoCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/KS2.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv")

control_epinano_coverage <- c1_epinano_coverage
control_epinano_coverage$c2_num_sites <- c2_epinano_coverage$num_sites
control_epinano_coverage$average_num_sites <- rowMeans(control_epinano_coverage[,c('num_sites', 'c2_num_sites')], na.rm=TRUE)

ks_epinano_coverage <- ks1_epinano_coverage
ks_epinano_coverage$ks2_num_sites <- ks2_epinano_coverage$num_sites
ks_epinano_coverage$average_num_sites <- rowMeans(ks_epinano_coverage[,c('num_sites', 'ks2_num_sites')], na.rm=TRUE)


# GRAPH
plot(
  type = "l",
  control_epinano_coverage$percent_through_genome,
  control_epinano_coverage$average_num_sites,
  xlab="% through genome",
  ylab="predicted number of m6a sites",
  col="blue"
)
lines(ks_epinano_coverage$percent_through_genome, ks_epinano_coverage$average_num_sites, col="red")

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}

# STATISTICAL TESTS
# plotting c1 against c2 num predicted sites should give a linear relationship. We can use a chow
# test to prove similarity between two series, and confirm dissimilartiy between knock sideways
#library(strucchange)

#sctest(control_epinano_coverage$num_sites ~ control_epinano_coverage$c2_num_sites, type="Chow")
#sctest(ks_epinano_coverage$num_sites ~ ks_epinano_coverage$ks2_num_sites, type="Chow")

# forget about the chow test, just plot the two and do a linear regression
#plot(control_epinano_coverage$num_sites, control_epinano_coverage$c2_num_sites)

scalar <- function(x) {
  x / sqrt(sum(x^2))
}

# and then to test against RRACH sites maybe normalise the two datasets then plot a linear regression
#plot(scalar(control_epinano_coverage$num_sites), scalar(rrach_genome_site_coverage$num_sites))
