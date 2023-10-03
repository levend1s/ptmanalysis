#!/usr/bin/env Rscript

# map epinano predicted sites to genes
library(stringr)
library(optparse)

# takes in a file of m6a site predictions and a gff file
# outputs a table of site predictions mapped to their features
# OPTIONAL takes in a limma-voom csv file and maps that to the genes aswell

args = commandArgs (trailingOnly=TRUE)
if (length(args) == 0 ) {
  stop ("-h, --help for usage")
}

option_list <- list (
  make_option (c("-p","--predictions_path"), type="character",
               help="path to file containing m6a predictions"),
  
  make_option (c("-t","--type"), type="character",
               help="type of predictions file (m6anet, epinano)"),
  
  make_option (c("-g","--gff_path"), type="character", 
               help="path to gff file containing features"),
  
  make_option (c("-d","--de_path"), type="character", 
               help="path to de file containing features (optional)"),
  
  make_option (c("-o","--out_path"), type="character", 
               help="path to write output to"),
  
  make_option (c("-v", '--verbose'), type="logical", default=0, action="store_true",
               help ="turn on verbose mode")
)

parser <- parse_args (OptionParser (option_list=option_list, usage="
  SitesToFeatures.R 0.1: Maps sites from prediction files to gff files. Can optionally map data from DE analysis (limma-voom csv)
  Currently supported types: m6anet, epinano"
))

if (!is.na(parser$predictions_path)) { 
  predictions_path <- parser$predictions_path
} else { 
  stop ('please provide the predictions file' ) 
}

if (!is.na(parser$type)) { 
  type <- parser$type
} else { 
  stop ('please provide the type of predictions file' ) 
}

if (!is.na(parser$de_path)) { 
  de_path <- parser$de_path
} else {
  de_path <- NA
}

if (!is.na(parser$gff_path)) { 
  gff_path <- parser$gff_path
} else { 
  stop ('please provide the gff file' ) 
}

if (!is.na(parser$out_path)) { 
  out_path <- parser$out_path
} else { 
  stop ('please provide the output file path' ) 
}

if (!is.na(parser$verbose)) { 
  verbose <- parser$verbose
} else { 
  verbose <- FALSE
}


# load data. NOTE these have header=FALSE because they have already been filtered to only include modified A's.
if (verbose) {
  write("loading data...", stdout())
}

if (type == "epinano") {
  predictions_data <- read.csv(predictions_path, header=FALSE, stringsAsFactors = TRUE)
  colnames(predictions_data) <- c("chr_pos", "ko_feature", "wt_feature", "delta_sum_err", "z_scores", "z_score_prediction")
  predictions_data[c("contig", "pos", "base", "strand")] <- str_split_fixed(predictions_data$chr_pos, " ", 4)
  predictions_data$pos <- as.numeric(predictions_data$pos) 
}
if (type == "m6anet") {
  predictions_data <- read.csv(predictions_path, header=TRUE, stringsAsFactors = TRUE)
}

if (!is.na(de_path)) {
  de <- read.csv(de_path, header=TRUE, stringsAsFactors = TRUE)
} else {
  de <- NA
}

gff_data <- read.delim(gff_path, header=FALSE, comment.char="#")
colnames(gff_data) <- c("contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# types that have "ID,description"
meta_feature_types <- list("protein_coding_gene", "pseudogene", "ncRNA_gene")
# types that have "ID,parent,..."
feature_types <- list("mRNA", "exon", "CDS", "pseudogenic_transcript", "five_prime_UTR", "three_prime_UTR", "ncRNA", "rRNA", "snoRNA", "tRNA", "snRNA")


# map each prediction site based off its contig, position and strand
# iterate through each m6a prediction site and find which genes it is in
mapToGFF <- function(data, gff, type) {
  mapped_m6as_list = list()
  
  for(i in 1:nrow(data)) {
    row <- data[i,]
    
    # get all pfal genes that the m6a site is in
    
    if (type == "epinano") {
      matches <- gff[
        gff$start <= row$pos 
        & gff$end >= row$pos
        & gff$contig == row$contig
        & gff$strand == row$strand,] 
      
      if (nrow(matches) == 0) {
        matches[1,] <- NA
      }
      
      df <- matches
      df$base_position <- row$pos
      df$strand <- row$strand
      df$contig <- row$contig
      df$z_scores <- row$z_scores
    }
    if (type == "m6anet") {
      matches <- gff[
        gff$start <= row$transcript_pos 
        & gff$end >= row$transcript_pos
        & gff$contig == row$transcript_id,] 
      
      if (nrow(matches) == 0) {
        matches[1,] <- NA
      }
      
      df <- matches
      df$base_position <- row$transcript_pos
      df$strand <- row$transcript_pos
      df$contig <- row$transcript_id
      df$probability_modified <- row$probability_modified
    }
    
    mapped_m6as_list[[i]] <- df
  }
  
  # now have m6a sites mapped to genes control
  mapped = do.call(rbind, mapped_m6as_list)
  
  return(mapped)
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

mapToDE <- function(prediction_table, de_table) {
  # the DE only looked at what featureCounts spits out, ie meta_features
  meta_features <- prediction_table[is.element(prediction_table$type, meta_feature_types),]
  prediction_table[c("logFC", "P.Value")] <- NA
  
  # iterate through each m6a prediction site, and add the FC and p value for the gene related to that site
  for(i in 1:nrow(de_table)) {
    row_de <- de_table[i,]

    if (nrow(prediction_table[!is.na(prediction_table$gene_id) & prediction_table$gene_id == row_de$gene_id,]) > 0) {
      prediction_table[!is.na(prediction_table$gene_id) & prediction_table$gene_id == row_de$gene_id,][c("logFC", "P.Value")] <- row_de[c("logFC", "P.Value")]
    }
  }
  
  return(prediction_table)
}

mapToDE <- function(prediction_table, de_table) {
  # the DE only looked at what featureCounts spits out, ie meta_features
  meta_features <- prediction_table[is.element(prediction_table$type, meta_feature_types),]
  prediction_table[c("logFC", "P.Value")] <- NA
  
  # iterate through each m6a prediction site, and add the FC and p value for the gene related to that site
  for(i in 1:nrow(de_table)) {
    row_de <- de_table[i,]
    
    if (nrow(prediction_table[!is.na(prediction_table$gene_id) & prediction_table$gene_id == row_de$gene_id,]) > 0) {
      prediction_table[!is.na(prediction_table$gene_id) & prediction_table$gene_id == row_de$gene_id,][c("logFC", "P.Value")] <- row_de[c("logFC", "P.Value")]
    }
  }
  
  return(prediction_table)
}

# split gff attributes
if (verbose) {
  write("splitting gff attributes...", stdout())
}
gff_data <- splitAttributes(gff_data)

# map predictions to GFF
if (verbose) {
  write("mapping predictions to gff...", stdout())
}

mapped_predictions <- mapToGFF(predictions_data, gff_data, type)
fraction_unmapped <- nrow(mapped_predictions[is.na(mapped_predictions),]) / nrow(mapped_predictions)
# drop all entries which have not been mapped to gene
mapped_predictions <- mapped_predictions[!is.na(mapped_predictions$gene_id),]

if (nrow(de) != 0) {
  # map predictions to DE analysis
  mapped_predictions <- mapToDE(mapped_predictions, de) 
}

if (verbose) {
  write("writing output...", stdout())
}
write.csv(mapped_predictions, out_path, row.names = FALSE)

if (verbose) {
  write("DONE.", stdout())
}
