#!/usr/bin/Rscript

library(AnnotationHub)
library(tidyverse)
library(regioneR)
library(ChIPpeakAnno)

args <- commandArgs(T)
species <- args[1]
dataProvider <- args[2]
version <- args[3]
bed_file <- args[4]
output_path <- args[5]

r_lib_paths <- .libPaths()
cache_path <- file.path(r_lib_paths[1], "AnnotationHub")
setAnnotationHubOption("CACHE", cache_path)
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c(species, dataProvider, version))
ahEdb <- ahDb[[names(ahDb)]]
trns <- transcripts(ahEdb)
	
##Get common biotypes
biotype_counts <- table(trns$tx_biotype)
biotype_counts_df <- data.frame(biotype_counts)
common_biotype <- biotype_counts_df[which(biotype_counts_df$Freq >= 1000),]
common_biotype_ar <- array(common_biotype$Var1)
	
variants <- toGRanges(bed_file)
annotation_result <- read_tsv(bed_file, col_names = c("chrom", "start", "end", "variant_id"))
	
##Calculate distance to different classes of genes
for(biotype in common_biotype_ar){
    temp <- trns[which(trns$tx_biotype == biotype),]
    annotatePeak = annotatePeakInBatch(variants, AnnotationData = temp, FeatureLocForDistance= "TSS", select = "arbitrary")
    name <- paste(biotype, "TSS", sep = "_")
    distance <- paste("distance_to", biotype, sep = "_")
    annotation_result[, name] <- annotatePeak$feature
    annotation_result[, distance] <- annotatePeak$distancetoFeature
}

write.table(annotation_result, output_path, row.names = FALSE, quote = FALSE, sep = "\t")


