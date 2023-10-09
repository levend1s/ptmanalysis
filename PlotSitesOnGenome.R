# read in gff file to calculate total size of genome bp
# assign resolution variable, ie window size to take avergae of, default should be genome size / 100 for 1% resolution
# count the number of m6a sites in each window size
# plot number of m6a sites against genome position

# --------------------- CHROMOSOME SUMMARY --------------------- 

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

gff <- read.delim(pfal_gff_path, header=FALSE, comment.char="#")
colnames(gff) <- c("contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff <- splitAttributes(gff)

# --------------------- ADD GENOMIC POSITION --------------------- 

addGenomicPosition <- function(d, type) {
  d <- d
  d$genomic_position <- NA
  
  if (type == "m6anet") {
    d$transcript_start_position <- NA
    d$contig <- NA
    d <- transform(d, transcript_start_position=gff[match(transcript_id, gff$ID), ]$start)
    # NOTE column transcript_position is the start position of the DRACH motiff. it is indexed +1, so transcript_position=37398 maps to contig position=37397
    # with this taken into account, the position of the A in the DRACH motif is only +1 above the column transcript_position
    d$position <- d$transcript_position+d$transcript_start_position+1
    d <- transform(d, contig=gff[match(transcript_id, gff$ID), ]$contig)
    d <- transform(d, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "plasmodb") {
    d <- transform(d, genomic_position=chromosomes[match(Genomic.Sequence.ID, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "epinanosvm") {
    d[c("start", "end")] <- str_split_fixed(d$Window, "-", 2)
    d$start <- as.numeric(d$start)
    d$end <- as.numeric(d$end)
    # always PLUS 2 from start, epinano doesnt show windows in reverse for - strand like plasmodb
    d$position <- d$start+2
    
    d <- transform(d, genomic_position=chromosomes[match(Ref, chromosomes$contig), ]$genomic_start)
  }
  else if(type == "epinano") {
    d <- transform(d, genomic_position=chromosomes[match(contig, chromosomes$contig), ]$genomic_start)
  }

  # add contig position to genomic position and subtract 1. This handles the case where genomic start = 1 and pos = 1 -> 1+1 - 1.
  d$genomic_position <- d$genomic_position + d$pos - 1

  return(d)
}

# --------------------- CALCULATE COVERAGE --------------------- 
# count how many DRACH sites there are every sliding_window_size bp genomic_position
num_data_points <- 100
sliding_window_size <- genome_size / num_data_points

calculateCoverage <- function(d, start, end) {
  cur_sliding_window_start <- start
  
  frac_num_points_to_perc <- num_data_points / 100
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
drach_positions_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/DRACH-DynSpansByMotifSearch_Summary.txt"
plasmodb_drachs <- read.table(file = drach_positions_path, sep = '\t', header = TRUE)
plasmodb_drachs$pos <- NA
plasmodb_drachs[plasmodb_drachs$Strand == "-",]$pos <- plasmodb_drachs[plasmodb_drachs$Strand == "-",]$Start-2
plasmodb_drachs[plasmodb_drachs$Strand == "+",]$pos <- plasmodb_drachs[plasmodb_drachs$Strand == "+",]$Start+2

# add genomic position to plasmodb_drachs
plasmodb_drachs <- addGenomicPosition(plasmodb_drachs, "plasmodb")
drach_genome_site_coverage <- calculateCoverage(plasmodb_drachs, 1, genome_size)

# ---------------------- data from DNA motif search using [GA][GA]AC[CAT] as pattern in Plasmodium falciparum 3D7
rrach_positions_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/RRACH-DynSpansByMotifSearch_Summary.txt"
plasmodb_rrachs <- read.table(file = rrach_positions_path, sep = '\t', header = TRUE)
plasmodb_rrachs$pos <- NA
plasmodb_rrachs[plasmodb_rrachs$Strand == "-",]$pos <- plasmodb_rrachs[plasmodb_rrachs$Strand == "-",]$Start-2
plasmodb_rrachs[plasmodb_rrachs$Strand == "+",]$pos <- plasmodb_rrachs[plasmodb_rrachs$Strand == "+",]$Start+2

# add genomic position to plasmodb_rrachs
plasmodb_rrachs <- addGenomicPosition(plasmodb_rrachs, "plasmodb")
rrach_genome_site_coverage <- calculateCoverage(plasmodb_rrachs, 1, genome_size)

# ---------------------- M6ANET ----------------------
generateM6AnetCoverage <- function(m6anet_data, start, end) {
  m6anet_probabilty_modified_threshold <- 0.8
  m6anet_data_filtered <- m6anet_data[m6anet_data$probability_modified > m6anet_probabilty_modified_threshold,]
  m6anet_genome_site_coverage <- calculateCoverage(m6anet_data_filtered, start, end)
  
  return(m6anet_genome_site_coverage)
}

c2_m6anet_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/inference_transcriptome/C2/data.site_proba.csv"
ks2_m6anet_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/inference_transcriptome/KS2/data.site_proba.csv"
c2_m6anet_data <- read.csv(c2_m6anet_path, header=TRUE, stringsAsFactors = TRUE)
ks2_m6anet_data <- read.csv(ks2_m6anet_path, header=TRUE, stringsAsFactors = TRUE)

c2_m6anet_data <- addGenomicPosition(c2_m6anet_data, "m6anet")
ks2_m6anet_data <- addGenomicPosition(ks2_m6anet_data, "m6anet")

c2_m6anet_coverage <- generateM6AnetCoverage(c2_m6anet_data, 1, genome_size)
ks2_m6anet_coverage <- generateM6AnetCoverage(ks2_m6anet_data, 1, genome_size)

control_m6anet_coverage <- c2_m6anet_coverage
#control_m6anet_coverage$c2_num_sites <- c2_m6anet_coverage$num_sites
#control_m6anet_coverage$average_num_sites <- rowMeans(control_m6anet_coverage[,c('num_sites', 'c2_num_sites')], na.rm=TRUE)

ks_m6anet_coverage <- ks2_m6anet_coverage
#ks_m6anet_coverage$ks2_num_sites <- ks2_m6anet_coverage$num_sites
#ks_m6anet_coverage$average_num_sites <- rowMeans(ks_m6anet_coverage[,c('num_sites', 'ks2_num_sites')], na.rm=TRUE)

# ---------------------- EPINANO ----------------------
generateEpinanoCoverage <- function(epinano_data, start, end) {
  epinano_data <- epinano_data[epinano_data$prediction == "mod",]
  epinano_genome_site_coverage <- calculateCoverage(epinano_data, start, end)
  
  return(epinano_genome_site_coverage)
}

c1_epinano_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/C1.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
c2_epinano_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/C2.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
ks1_epinano_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/KS1.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"
ks2_epinano_path <- "/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/epinano/PredictionSVM/KS2.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"

c1_epinano_data <- read.csv(c1_epinano_path, header=TRUE, stringsAsFactors = TRUE)
c2_epinano_data <- read.csv(c2_epinano_path, header=TRUE, stringsAsFactors = TRUE)
ks1_epinano_data <- read.csv(ks1_epinano_path, header=TRUE, stringsAsFactors = TRUE)
ks2_epinano_data <- read.csv(ks2_epinano_path, header=TRUE, stringsAsFactors = TRUE)

c1_epinano_data <- addGenomicPosition(c1_epinano_data, "epinanosvm")
c2_epinano_data <- addGenomicPosition(c2_epinano_data, "epinanosvm")
ks1_epinano_data <- addGenomicPosition(ks1_epinano_data, "epinanosvm")
ks2_epinano_data <- addGenomicPosition(ks2_epinano_data, "epinanosvm")

c1_epinano_coverage <- generateEpinanoCoverage(c1_epinano_data, 1, genome_size)
c2_epinano_coverage <- generateEpinanoCoverage(c2_epinano_data, 1, genome_size)
ks1_epinano_coverage <- generateEpinanoCoverage(ks1_epinano_data, 1, genome_size)
ks2_epinano_coverage <- generateEpinanoCoverage(ks2_epinano_data, 1, genome_size)

control_epinano_coverage <- c1_epinano_coverage
control_epinano_coverage$c2_num_sites <- c2_epinano_coverage$num_sites
control_epinano_coverage$average_num_sites <- rowMeans(control_epinano_coverage[,c('num_sites', 'c2_num_sites')], na.rm=TRUE)

ks_epinano_coverage <- ks1_epinano_coverage
ks_epinano_coverage$ks2_num_sites <- ks2_epinano_coverage$num_sites
ks_epinano_coverage$average_num_sites <- rowMeans(ks_epinano_coverage[,c('num_sites', 'ks2_num_sites')], na.rm=TRUE)

# ---------------------- GRAPH ---------------------- 
scalar <- function(x) {
  x / sqrt(sum(x^2))
}

# POTENTIAL SITES
plot(
  type = "l",
  drach_genome_site_coverage$percent_through_genome,
  drach_genome_site_coverage$num_sites,
  xlab="% through genome",
  ylab="number of DRACH/RRACH sites",
  ylim=c(0, 8500)
)

lines(rrach_genome_site_coverage$percent_through_genome, rrach_genome_site_coverage$num_sites, col="red")

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}

plot(
  type = "l",
  control_epinano_coverage$percent_through_genome,
  scalar(control_epinano_coverage$average_num_sites),
  xlab="% through genome",
  ylab="number predicted sites in control across tools",
  ylim=c(0,1),
  col="blue"
)
lines(control_m6anet_coverage$percent_through_genome, scalar(control_m6anet_coverage$num_sites), col="red")

for(i in 1:nrow(chromosomes)) {
  abline(v=chromosomes[i,]$genomic_start_percentage, col="gray", lty=2)
}

# do it just for chomosome x
chromosome_x <- 14
chr1_c2_m6anet_coverage <- generateM6AnetCoverage(c2_m6anet_data, chromosomes[chromosome_x,]$genomic_start, chromosomes[chromosome_x,]$genomic_start + chromosomes[chromosome_x,]$end)
#ks1_m6anet_coverage <- generateM6AnetCoverage("/Users/joshualevendis/Desktop/Biomedical\ Research\ Project/m6anet/ks_1_m6anet_output_inference/ks1_data.site_proba.csv")
chr1_ks2_m6anet_coverage <- generateM6AnetCoverage(ks2_m6anet_data, chromosomes[chromosome_x,]$genomic_start, chromosomes[chromosome_x,]$genomic_start + chromosomes[chromosome_x,]$end)

# GRAPH
plot(
  type = "l",
  chr1_c2_m6anet_coverage$percent_through_genome,
  chr1_c2_m6anet_coverage$num_sites,
  xlab="% through genome",
  ylab="number of DRACH sites",
  col="blue"
)
lines(chr1_ks2_m6anet_coverage$percent_through_genome, chr1_ks2_m6anet_coverage$num_sites, col="red")

# STATISTICAL TESTS
# plotting c1 against c2 num predicted sites should give a linear relationship. We can use a chow
# test to prove similarity between two series, and confirm dissimilartiy between knock sideways
#library(strucchange)

#sctest(control_epinano_coverage$num_sites ~ control_epinano_coverage$c2_num_sites, type="Chow")
#sctest(ks_epinano_coverage$num_sites ~ ks_epinano_coverage$ks2_num_sites, type="Chow")

# forget about the chow test, just plot the two and do a linear regression
#plot(control_epinano_coverage$num_sites, control_epinano_coverage$c2_num_sites)

# and then to test against RRACH sites maybe normalise the two datasets then plot a linear regression
#plot(scalar(control_epinano_coverage$num_sites), scalar(rrach_genome_site_coverage$num_sites))
