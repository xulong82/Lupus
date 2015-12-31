library(affy)
library(oligo)
library(limma)
library(DESeq2)
library(pheatmap)
library(VennDiagram)
library(gplots)
library(dplyr)
library(xlsx)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

yaa <- list()

# BXSB_YAA: RNA-seq
eset <- read.delim("Greg/counts.txt", stringsAsFactors = F)
eset <- data.frame(row.names = eset$Gene, eset[grep("BXSB", colnames(eset))])
eset <- eset[apply(eset, 1, function(x) max(x) > 10 & sum(x > 0) > 2), ] 

group <- factor(gsub(".*(B6|Yaa).*", "\\1", colnames(eset)), levels = c("B6", "Yaa"))
colData <- data.frame(row.names = colnames(eset), condition = group)
dds <- DESeqDataSetFromMatrix(countData = eset, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
plotMA(res, main = "DESeq2", ylim = c(-2,2))
summary(res)

resSig <- subset(res, padj < 0.05) %>% as.data.frame
resSig <- subset(res, padj < 0.1) %>% as.data.frame
resSig <- subset(res, pvalue < 0.05) %>% as.data.frame  # Saved
yaa$bxsb$sig <- resSig
yaa$bxsb$up <- up <- rownames(resSig)[resSig$log2FoldChange > 0]
yaa$bxsb$down <- down <- rownames(resSig)[resSig$log2FoldChange < 0]

# B6_YAA: MicroArray
cel.files <- read.celfiles(list.celfiles("bRNA/B6/CEL", full.name = T))
norm.rma <- rma(cel.files, target = "core")  # RMA, gene level summarization
show(norm.rma)
eset <- exprs(norm.rma)
colnames(eset) <- gsub("GC_Gene1-0ST_(.*)_430.*", "\\1", colnames(eset))

sampleInf <- read.xlsx("bRNA/B6/Sample.xlsx", sheetName = "Sheet1", stringsAsFactors = F)
sampleInf <- sampleInf[match(colnames(eset), sampleInf$Sample.Name), ]
sampleInf$Genotype[grep("wild type", sampleInf$Genotype)] <- "WT"
sampleInf$Genotype[grep("other", sampleInf$Genotype)] <- "MT"

mogene <- read.delim("bRNA/B6/mogene2symbol.map", sep = "\t", stringsAsFactors = F)
eset <- eset[rownames(eset) %in% mogene[, 1], ]
symbol <- mogene$Associated.Gene.Name[match(rownames(eset), mogene$Affy.MoGene.probeset)]
eset <- apply(eset, 2, function(x) tapply(x, symbol, max))
eset <- eset[apply(eset, 1, function(x) min(x) > log2(120)), ]

fit <- apply(eset, 1, function(x) lm(x ~ sampleInf$Genotype))
fit.r2 <- sapply(fit, function (x) summary(x)$r.squared)
fit.fs <- sapply(fit, function (x) summary(x)$fstatistic)
fit.et <- sapply(fit, function (x) summary(x)$coefficients[2, 1])
fit.pv <- apply(fit.fs, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F)) 
fit.qv <- qvalue(fit.pv)

yaa$b6$sig <- fit.sig <- fit[fit.pv < 0.05]
yaa$b6$up <- up <- names(fit)[fit.pv < 0.05 & fit.et > 0]
yaa$b6$down <- down <- names(fit)[fit.pv < 0.05 & fit.et < 0]

# ---
load("data/yaaList.rdt")
for(obj in names(yaa)) assign(obj, yaa[[obj]])

venn.diagram(list(B6_up = b6$down, B6_down = b6$up, BXSB_up = bxsb$down, BXSB_down = bxsb$up), 
             file = "Figure/venn.tiff", fill = c("firebrick1", "darkorchid2", "chartreuse3", "dodgerblue3"))

bxsb1 <- setdiff(bxsb$down, c(b6$up, b6$down))
bxsb2 <- intersect(bxsb$up, c(b6$up, b6$down))

source("../X/function.R")
source("../SCR/kegg.R")
gk1 <- kegg(bxsb1)
gk2 <- kegg(bxsb2)

keggFind("pathway", "Lupus")  # KEGG: SLE
write(keggGet("mmu05322", "kgml"), file = "mmu05322.xml")  # Cytoscape
writePNG(keggGet("mmu05322", "image"), "mmu05322.png")
writePNG(keggGet("hsa05322", "image"), "hsa05322.png")

# IL21: RNA-seq
load("data/profile12.rdt")  # profile1

# --- INTEGRATE GENE EXPRESSION EXPERIMENTS, MESSY CODE, CAN BE DELETED

# Structure of the return, y (A better algorithm?)
#      NN  NP  PP
# -------------------
#  B6  1   2   3
#  YAA 4   5   6

myCluster <- function(x) {
  y <- rep(0, nrow(x))
  for (i in 1:nrow(x)) {
    m1 <- max(x[i, 1:3])
    m2 <- max(x[i, 4:5])
    if (m1 == x[i, 1]) {
      y[i] <- ifelse(m2 == x[i, 4], 1, 4)
    } else if (m1 == x[i, 2]) {
      y[i] <- ifelse(m2 == x[i, 4], 2, 5)
    } else {
      y[i] <- ifelse(m2 == x[i, 4], 3, 6)
    }
  }
  return(y)
}

# Significant genes in each experiment
list.bx <- sig.deseq[sig.deseq$id %in% rownames(il21), ]$id  # Sig genes in BXSB strain
list.b6 <- rownames(probeset.list)[rownames(probeset.list) %in% rownames(il21)]  # Sig genes in B6
list.il <- rownames(il21.2)

# Filter: high genes in each experiment 
hg.il <- rownames(il21[apply(il21, 1, max) > 100, ])  # High genes in IL21
hg.bx <- rownames(bxsb.counts.m[apply(bxsb.counts.m, 1, max) > 100, ])  # High genes in BXSB
hg.b6 <- rownames(b6.yaa.m[apply(b6.yaa.m, 1, max) > 5, ])  # High genes in B6
list.bx <- intersect(list.bx, hg.bx)
list.b6 <- intersect(list.b6, hg.b6)
list.il <- intersect(list.il, hg.il)

# Correlation on significant genes 
data.b6 <- b6.yaa.m[rownames(b6.yaa.m) %in% list.b6, ]
data.b6 <- data.b6[order(rownames(data.b6)), ]
data.bx <- bxsb.counts.m[rownames(bxsb.counts.m) %in% list.bx, ]
data.bx <- data.bx[order(rownames(data.bx)), ]

data.il <- il21[rownames(il21) %in% list.b6, ]
data.il <- il21[rownames(il21) %in% list.bx, ]
data.il <- data.il[order(rownames(data.il)), ]

data.cb <- cbind(data.b6, data.il)
data.cb <- cbind(data.bx, data.il)
pheatmap(cor(data.cb, method = "spearman"), fontsize = 16, fontsize_number = 10, display_number = T, 
         treeheight_row = 150, treeheight_col = 0, cellwidth = 38, cellheight = 40)

bx.il <- intersect(list.bx, list.il)
il21.3 <- il21.2[rownames(il21.2) %in% bx.il, ]
il21.3 <- il21.3[order(rownames(il21.3)), ][, 6:8]
bxsb.2 <- sig.deseq[sig.deseq$id %in% bx.il, ]
bxsb.2 <- bxsb.2[order(bxsb.2$id), ][, 3:4]
data.bx.il <- cbind(il21.3, bxsb.2)
colnames(data.bx.il) <- c("VNIN", "VNIP", "VPIP", "B6", "YAA")

b6.il <- intersect(list.b6, list.il)
il21.4 <- il21.2[rownames(il21.2) %in% b6.il, ]
il21.4 <- il21.4[order(rownames(il21.4)), ][, 6:8]
b6.2 <- b6.yaa.m[rownames(b6.yaa.m) %in% b6.il, ]
b6.2 <- b6.2[order(rownames(b6.2)), ]
data.b6.il <- cbind(il21.4, b6.2)
colnames(data.b6.il) <- c("VNIN", "VNIP", "VPIP", "B6", "YAA")

clus.bx.il <- myCluster(data.bx.il)  # up/down clusters
clus.b6.il <- myCluster(data.b6.il)

table(clus.bx.il)
table(clus.b6.il)

par(mfrow = c(1, 2))
venn(list(BXSB = list.bx, IL21 = list.il))
venn(list(B6 = list.b6, IL21 = list.il))
par(mfrow = c(1, 1))
venn.diagram(list(B6 = list.b6, IL21 = list.il, BXSB = list.bx), 
             file = "1.tiff", fill = c("red", "green", "blue"))

table.bx <- matrix(c(200, 46, 108, 211, 70, 245), nrow = 2)
table.b6 <- matrix(c(27, 32, 9, 75, 13, 78), nrow = 2)
fisher.test(table.bx)
fisher.test(table.b6)
# chisq.test(table.b6[1, ])
# chisq.test(table.b6[2, ])

bx.il.id <- matrix(c(rownames(data.bx.il)[clus.bx.il == 6], rep("NULL", 1)), ncol = 10)
dimnames(bx.il.id) <- list(c(1:18), c(1:10))
b6.il.id <- matrix(c(rownames(data.b6.il)[clus.b6.il == 6], rep("NULL", 7)), ncol = 10)
dimnames(b6.il.id) <- list(c(1:8), c(1:10))
textplot(bx.il.id, cex = .7)
textplot(b6.il.id, cex = .7)

bx.b6.id <- intersect(bx.il.id, b6.il.id)
bx.b6.id <- matrix(c(bx.b6.id, rep("NULL", 0)), ncol = 7)
dimnames(bx.b6.id) <- list(c(1:5), c(1:7))
venn.diagram(list(BXSB = 1:179, B6 = 1:73 + (179 - 35)), file = "3.tiff", fill = c("red", "green"))
textplot(bx.b6.id)

yaa <- c("Msl3", "Tlr7", "Tmsb4x", "Rab9", "Prps2", "Trappc2", "Mid1", "Arhgap6")
intersect(list.b6, yaa)
intersect(list.bx, yaa)
intersect(list.il, yaa)

tfdb <- read.delim("~/Dropbox/Lupus/TFdb.Riken.txt", header = F)  # 1675 TFs from Riken TFdb
intersect(tfdb$V1, rownames(data.bx.il)[clus.bx.il == 6])

gwas <- read.delim("~/Dropbox/Lupus/gwas.genes2.txt", header = F)  # MGI and GWAS hints of SLE

intersect(bx.b6.id, tfdb$V1)
intersect(bx.b6.id, gwas$V1)

intersect(list.bx, gwas$V1)
intersect(list.b6, gwas$V1)
intersect(list.il, gwas$V1)
gwas.list <- c("Rasgrp1", "Stat6", "Il10", "Trove2")
gwas.il <- il21[rownames(il21) %in% gwas.list, ]
gwas.bx <- bxsb.counts.m[rownames(bxsb.counts.m) %in% gwas.list, ]
gwas.b6 <- b6.yaa.m[rownames(b6.yaa.m) %in% gwas.list, ]

plot(rowMeans(b6.yaa[, c(1:3)]), rowMeans(b6.yaa[, c(4:5)]), xlab = "WT", ylab = "Yaa")
pheatmap(cor(b6.yaa), fontsize = 16, fontsize_number = 10, display_number = T, 
         treeheight_row = 150, treeheight_col = 0, cellwidth = 38, cellheight = 40)
