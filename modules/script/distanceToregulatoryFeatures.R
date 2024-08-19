#!/usr/bin/Rscript

library(tidyverse)
library(regioneR)
library(ChIPpeakAnno)

args <- commandArgs(T)
bed_file <- args[1]
output_path <- args[2]
enhancer_path <- args[3]
TF_path <- args[4]
CTCF_path <- args[5]
promoter_path <- args[6]

annotation_result <- read_tsv(bed_file, col_names = c("chrom", "start", "end", "variant_id"))
variants <- toGRanges(bed_file)

enhancer_gr <- toGRanges(enhancer_path)
TF_gr <- toGRanges(TF_path)
CTCF_gr <- toGRanges(CTCF_path)
promoter_gr <- toGRanges(promoter_path)

annotatePeak_enhancer = annotatePeakInBatch(variants, AnnotationData = enhancer_gr, FeatureLocForDistance= "TSS", select = "arbitrary")
annotatePeak_TF = annotatePeakInBatch(variants, AnnotationData = TF_gr, FeatureLocForDistance= "TSS", select = "arbitrary")
annotatePeak_CTCF = annotatePeakInBatch(variants, AnnotationData = CTCF_gr, FeatureLocForDistance= "TSS", select = "arbitrary")
annotatePeak_promoter = annotatePeakInBatch(variants, AnnotationData = promoter_gr, FeatureLocForDistance= "TSS", select = "arbitrary")

annotation_result[, "distance_to_enhancers"] <- annotatePeak_enhancer$distancetoFeature
annotation_result[, "distance_to_TF"] <- annotatePeak_TF$distancetoFeature
annotation_result[, "distance_to_CTCF"] <- annotatePeak_CTCF$distancetoFeature
annotation_result[, "distance_to_promoter"] <- annotatePeak_promoter$distancetoFeature
write.table(annotation_result, output_path, row.names = FALSE, quote = FALSE, sep = "\t")

