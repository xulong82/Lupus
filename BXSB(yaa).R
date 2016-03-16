library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Gviz)

options(stringsAsFactors = F)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

# yaa and par loci

axTrack <- GenomeAxisTrack()
gTrackX <- IdeogramTrack(genome = "mm10", chromosome = "chrX")
gTrackY <- IdeogramTrack(genome = "mm10", chromosome = "chrY")

par <- list(chrX = c(169969759, 170931299), chrY = c(90745845, 91644698))
yaa <- list(chrX = c(161031299, 171031299), chrY = c(81644698, 91644698)) # 10Mbp

refGenes <- function(chr, from, to) UcscTrack(genome="mm10", chromosome=chr, track="refGene", from=from, to=to,
  trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", name="Ensembl Genes", fill="red",
  gene="name", symbol="name2", transcript="name", strand = "strand")

pdf("./wes/par.pdf", width = 12, height = 4)

txTrack <- refGenes("chrX", par$chrX[1], par$chrX[2])
plotTracks(list(gTrackX, axTrack, txTrack), from=par$chrX[1], to=par$chrX[2], transcriptAnnotation="symbol", stackHeight=0.2)
txTrack <- refGenes("chrY", par$chrY[1], par$chrY[2])
plotTracks(list(gTrackY, axTrack, txTrack), from=par$chrY[1], to=par$chrY[2], transcriptAnnotation="symbol", stackHeight=0.2)
dev.off()

pdf("./wes/yaa.pdf", width = 12, height = 4)

txTrack <- refGenes("chrX", yaa$chrX[1], yaa$chrX[2])
hlTrack <- HighlightTrack(trackList = txTrack, start = par$chrX[1], end = par$chrX[2], chromosome = "chrX", fill = "Skyblue")
plotTracks(list(gTrackX, axTrack, hlTrack), from=yaa$chrX[1], to=yaa$chrX[2], transcriptAnnotation="symbol", stackHeight=0.2)
symbol(txTrack) %>% unique

txTrack <- refGenes("chrY", yaa$chrY[1], yaa$chrY[2])
hlTrack <- HighlightTrack(trackList = txTrack, start = par$chrY[1], end = par$chrY[2], chromosome = "chrY", fill = "Skyblue")
plotTracks(list(gTrackY, axTrack, hlTrack), from=yaa$chrY[1], to=yaa$chrY[2], transcriptAnnotation="symbol", stackHeight=0.2)

dev.off()

pdf("./wes/chrY.pdf", width = 12, height = 4)

txTrack <- refGenes("chrY", 1, yaa$chrY[2])
hlTrack <- HighlightTrack(trackList = list(txTrack), start = par$chrY[1], end = par$chrY[2], chromosome = "chrY")
plotTracks(list(gTrackY, axTrack, hlTrack), from=1, to=yaa$chrY[2], transcriptAnnotation="symbol", stackHeight=0.2)

dev.off()
