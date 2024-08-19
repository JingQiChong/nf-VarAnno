#!/usr/bin/Rscript

library(AnnotationHub)
library(tidyverse)
library(regioneR)

args <- commandArgs(T)
species <- args[1]
dataProvider <- args[2]
version <- args[3]
bed_file <- args[4]
output_path <- args[5]

## Get the gene density in tilewidth
r_lib_paths <- .libPaths()
cache_path <- file.path(r_lib_paths[1], "AnnotationHub")
setAnnotationHubOption("CACHE", cache_path)
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c(species, dataProvider, version))
ahEdb <- ahDb[[names(ahDb)]]
human.genes <- genes(ahEdb)
human.windows<-tileGenome(seqinfo(human.genes), tilewidth=1000000, cut.last.tile.in.chrom=T)
human.windows$numofGenes<-countOverlaps(human.windows, human.genes)

## variants_gr and human.windows have different seqlevels (name of chromosome), need same seqlevels
variants_gr <- toGRanges(bed_file)
newStyle <- mapSeqlevels(seqlevels(variants_gr), "NCBI")
variants_gr <- renameSeqlevels(variants_gr, newStyle)

## Get the gene desity of the each variant
overlaps <- findOverlaps(human.windows, variants_gr)
gene_density <- human.windows[queryHits(overlaps)]$numofGenes

annotation_result <- read_tsv(bed_file, col_names = c("chrom", "start", "end", "variant_id"))
annotation_result$gene_density <- gene_density
write.table(annotation_result, output_path, row.names = FALSE, quote = FALSE, sep = "\t")