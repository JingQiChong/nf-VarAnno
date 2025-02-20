#!/usr/bin/Rscript

library(AnnotationHub)
library(tidyverse)
library(regioneR)
library(ChIPpeakAnno)
library(BiocFileCache)

args <- commandArgs(T)
species <- args[1]
dataProvider <- args[2]
version <- args[3]
bed_file <- args[4]
output_path <- args[5]

# Update AnnotationHub cache path

moveFiles<-function(package){
    olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
    newdir <- tools::R_user_dir(package, which="cache")
    dir.create(path=newdir, recursive=TRUE)   
    files <- list.files(olddir, full.names =TRUE)
    moveres <- vapply(files,
    FUN=function(fl){
        filename = basename(fl)
        newname = file.path(newdir, filename)
        file.rename(fl, newname)
    },
    FUN.VALUE = logical(1))
    if(all(moveres)) unlink(olddir, recursive=TRUE)
}
package="AnnotationHub"
moveFiles(package)

ah <- AnnotationHub()
ahDb <- query(ah, pattern = c(species, dataProvider, version))
ahEdb <- ahDb[[names(ahDb), force=TRUE]]
trns <- transcripts(ahEdb)
	
# Get common biotypes
biotype_counts <- table(trns$tx_biotype)
biotype_counts_df <- data.frame(biotype_counts)
common_biotype <- biotype_counts_df[which(biotype_counts_df$Freq >= 1000),]
common_biotype_ar <- array(common_biotype$Var1)
	
variants <- toGRanges(bed_file)
annotation_result <- read_tsv(bed_file, col_names = c("chrom", "start", "end", "variant_id"))
	
# Calculate distance to different classes of genes
for(biotype in common_biotype_ar){
    temp <- trns[which(trns$tx_biotype == biotype),]
    annotatePeak = annotatePeakInBatch(variants, AnnotationData = temp, FeatureLocForDistance= "TSS", select = "arbitrary")
    name <- paste(biotype, "TSS", sep = "_")
    distance <- paste("distance_to", biotype, sep = "_")
    annotation_result[, name] <- annotatePeak$feature
    annotation_result[, distance] <- annotatePeak$distancetoFeature
}

write.table(annotation_result, output_path, row.names = FALSE, quote = FALSE, sep = "\t")


