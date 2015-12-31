library(dplyr)
library(ggplot2)
library(genefilter)
library(Gviz)
library(Biobase)
library(biomaRt)

options(stringsAsFactors = F)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

load("./BXSB(rnaseq)/gene_expression_tpm.rdt")
data <- gene_expression_tpm[rowMax(gene_expression_tpm) > 5, ] %>% as.data.frame
data <- data[! rowMax(as.matrix(data)) > 1e4, ]
data <- sweep(data, 2, colSums(data), "/") * 1e6 
data$BXSB_avg <- rowMeans(data[group == "BXSB"])
data$BXSB_B6_avg <- rowMeans(data[group == "BXSB_B6"])

Il21 <- list(chr = "chr3", tss = 37222759)
range <- list(chr = Il21$chr, start = Il21$tss - 1e7, end = Il21$tss + 1e7)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attribute <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")
filters <- c("chromosome_name", "start", "end")

genes <- getBM(attribute, filters, values = list(3, range$start, range$end), mart)
genes <- filter(genes, ensembl_gene_id %in% rownames(data))
genes <- cbind(genes, data[genes$ensembl_gene_id, c("BXSB_avg", "BXSB_B6_avg")])
genes$chromosome_name <- paste0("chr", genes$chromosome_name)
genes$external_gene_name

rnaGR <- makeGRangesFromDataFrame(genes, keep.extra.columns = T, 
  seqnames.field = "chromosome_name", start.field = "start_position", end.field = "end_position")

rnaTrack <- DataTrack(rnaGR, name = "RNA-seq", groups=rep(c("BXSB", "BXSB_B6"), each = 1))
plotTracks(rnaTrack, from = range$start, to = range$end, type = "p", legend = T, lwd = 2, cex.legend = 1)

var <- read.delim("./WES/BXSB_SNP_filtered.vcf", skip = 23) # filter required
names(var)[1] <- "CHROM"
var$CHROM <- factor(var$CHROM, levels = paste0("chr", c(1:19, "X", "Y")))
table(var$CHROM)
var_Il21 = filter(var, CHROM == range$chr & POS > range$start & POS < range$end)

varGR <- makeGRangesFromDataFrame(var_Il21, seqnames.field = "CHROM", start.field = "POS", end.field = "POS")
varTrack <- AnnotationTrack(varGR, name = "var", col = "red")
plotTracks(varTrack, from = range$start, to = range$end, col = "red")

geneGR <- with(genes, GRanges(chromosome_name, IRanges(start_position, end_position), id = external_gene_name))
grTrack <- AnnotationTrack(geneGR, name = "Genes", showFeatureId = T, fill = "darkgreen", fontcolor.feature = "red")
plotTracks(grTrack, stacking = "full")

pdf("pdf/gviz.pdf", width = 9, height = 5)
plotTracks(list(rnaTrack, varTrack), legend = T, lwd = 2)
dev.off()

plot(genes$start_position, genes$BXSB_avg)

ggplot(var_Il21, aes(x = POS)) + geom_bar(binwidth = 1e5, fill = "red") +
  theme_bw() + xlab("") + ylab("Count") + coord_cartesian(ylim=c(0, 9))
