#!/usr/bin/Rscript

library(AnnotationHub)
library(tidyverse)
library(regioneR)
library(BiocFileCache)

args <- commandArgs(T)
species <- args[1]
dataProvider <- args[2]
version <- args[3]
bed_file <- args[4]
output_path <- args[5]

# Update the AnnotationHub cache path
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
human.genes <- genes(ahEdb)
human.windows<-tileGenome(seqinfo(human.genes), tilewidth=1000000, cut.last.tile.in.chrom=T)
human.windows$numofGenes<-countOverlaps(human.windows, human.genes)

variants_gr <- toGRanges(bed_file)
newStyle <- mapSeqlevels(seqlevels(variants_gr), "NCBI")
variants_gr <- renameSeqlevels(variants_gr, newStyle)

overlaps <- findOverlaps(human.windows, variants_gr)
gene_density <- human.windows[queryHits(overlaps)]$numofGenes

annotation_result <- read_tsv(bed_file, col_names = c("chrom", "start", "end", "variant_id"))
annotation_result$gene_density <- gene_density
write.table(annotation_result, output_path, row.names = FALSE, quote = FALSE, sep = "\t")