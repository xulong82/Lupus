library(dplyr)
library(genefilter)
library(Biobase)
library(ggplot2)
library(biomaRt)
library(qvalue)
library(DESeq2)
library(preprocessCore)
library(quantro)
library(matrixStats)
library(Hmisc)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")
load("bRNA/gene_expression_tpm.rdt")
source("../../X/function.R")

data <- gene_expression_tpm[rowMax(gene_expression_tpm) > 10, ] %>% as.data.frame
group <- factor(names(data), levels = c("BXSB", "BXSB_B6"))

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attribute <- c("ensembl_gene_id", "external_gene_name", "description")

pdf("bRNA/Igkc.pdf", width = 8, height = 6)
par(mar = c(10, 4, 4, 2))
matboxplot(data, groupFactor = group)
dev.off()

outlier <- data[rowMax(as.matrix(data)) > 1e4, ] # ENSMUSG00000076609/Igkc is twice in BXSB
mart1 <- getBM(attribute, "ensembl_gene_id", rownames(outlier), mart)
mart1$description <- gsub("_\\[.*", "", mart1$description)

data <- data[! rowMax(as.matrix(data)) > 1e4, ]
data <- sweep(data, 2, colSums(data), "/") * 1e6 

pdf("bRNA/data.pdf", width = 8, height = 6)
par(mar = c(10, 4, 4, 2))
matboxplot(log2(data + 1), groupFactor = group)
dev.off()

objectMedians <- colMedians(as.matrix(data))
normalize.median <- sweep(data, 2, objectMedians, FUN = "-")
colMedians(as.matrix(normalize.median))
colSums(normalize.median)
normalize.quantile <- normalize.quantiles(as.matrix(data), copy = TRUE)
dimnames(normalize.quantile) <- dimnames(data)

par(mfrow = c(1, 2))
qqplot(data[, 1], data[, 7], xlab = "BXSB_B6", ylab = "BXSB")
abline(0, 1)
qqplot(normalize.median[, 1], normalize.median[, 7], xlab = "BXSB_B6", ylab = "BXSB")
abline(0, 1)

ttest <- rowttests(as.matrix(log2(data + 1)), group) # ttests

data$BXSB_avg <- rowMeans(data[group == "BXSB"])
data$BXSB_B6_avg <- rowMeans(data[group == "BXSB_B6"])
data$log2fold <- with(data, log2(BXSB_B6_avg + 1) - log2(BXSB_avg + 1))
data <- cbind(data, ttest)
data$q.value <- qvalue(data$p.value)$qvalues
data_select <- data[abs(data$log2fold) > 0.3 & data$q.value < 0.05, ]
ttest$select = "N"
ttest$select[rownames(ttest) %in% rownames(data_select)] = "S"
# !!! transcriptome of BXSB_B6 and BXSB so different!

pdf("bRNA/vocano.pdf", width = 6, height = 4)
ggplot(ttest, aes(x = dm, y = -log10(p.value))) + 
  geom_point(aes(color = as.factor(select)), size = 1) +
  scale_color_manual(values = c("grey30", "firebrick1")) +
  theme_bw() + xlab("Effect") + ylab("-log10(pvalue)") +
  theme(legend.position = "none")
dev.off()

mart2 <- getBM(attribute, "ensembl_gene_id", rownames(data_select), mart)
symbol <- mart2$external_gene_name %>% unique
gk <- myGK(symbol)
gk$KEGG[1:10, ]
sapply(gk$GO, function(x) x$Term[1:20])
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

RNA = list(data = data, symbol = symbol, gk = gk)
bxsbList = list(RNA = RNA)

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
group <- factor(sampleInf$Genotype, levels = c("WT", "MT"))

mogene <- read.delim("bRNA/B6/mogene2symbol.map", sep = "\t", stringsAsFactors = F)
eset <- eset[rownames(eset) %in% mogene[, 1], ]
symbol <- mogene$Associated.Gene.Name[match(rownames(eset), mogene$Affy.MoGene.probeset)]
eset <- apply(eset, 2, function(x) tapply(x, symbol, max))
eset <- eset[apply(eset, 1, function(x) min(x) > log2(120)), ] %>% as.data.frame

ttest <- rowttests(as.matrix(eset), group) # ttests

eset$WT_avg <- rowMeans(eset[group == "WT"])
eset$MT_avg <- rowMeans(eset[group == "MT"])
eset$log2fold <- with(eset, log2(MT_avg) - log2(WT_avg))
eset <- cbind(eset, ttest)
eset$q.value <- qvalue(eset$p.value)$qvalues

eset_select <- eset[eset$p.value < 0.05, ]
gk <- myGK(rownames(eset_select))
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])
bxsbList$B6Array <- eset

ttest$select = "N"
ttest$select[rownames(ttest) %in% rownames(eset_select)] = "S"
pdf("bRNA/vocano_b6.pdf", width = 6, height = 4)
ggplot(ttest, aes(x = dm, y = -log10(p.value))) + 
  geom_point(aes(color = as.factor(select)), size = 1) +
  scale_color_manual(values = c("grey30", "firebrick1")) +
  theme_bw() + xlab("Effect") + ylab("-log10(pvalue)") +
  theme(legend.position = "none")
dev.off()

# Signature genes
yaa_bxsb <- bxsbList$RNA$symbol
yaa_b6 <- rownames(eset_select) 
intersect(yaa_bxsb, yaa_b6)
setdiff(yaa_bxsb, yaa_b6)
geneList <- list(B6 = setdiff(yaa_b6, yaa_bxsb), BXSB = setdiff(yaa_bxsb, yaa_b6), YAA = intersect(yaa_bxsb, yaa_b6))
venn.diagram(geneList, imagetype = "png", file = "bRNA/venn.png", width = 2000, height = 2000)

gk = myGK(geneList$B6)
gk = myGK(geneList$BXSB)
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

load("./Immgen/immgenList.rdt") # module
module <- lapply(module, function(x) capitalize(tolower(x)))
bg <- unique(c(rownames(eset), yaa_bxsb, unlist(module)))

myhyper <- function(g1, g2) {  # Hypergeometric
  if(length(intersect(g1, g2)) == 0) return(1)
  1 - phyper(length(intersect(g1, g2)) - 1, length(g2), length(setdiff(bg, g2)), length(g1))
}  # Pr(count >= length(intersect(g1, g2)))

sapply(module, length)
(int = sapply(module, function(x) length(intersect(geneList$BXSB, x))))
(enr = sapply(module, function(x) myhyper(geneList$BXSB, x)))
df = as.data.frame(cbind(int, enr))
df$module1 = 1:81
df$module2 = rownames(df)
df$text = NULL
df$text[df$enr < 0.01] = which(df$enr < 0.01)

pdf("bRNA/immgen.pdf", width = 8, height = 5)
ggplot(df, aes(x = module1, y = -log10(enr), label = text)) +  
  geom_point(aes(color = module2, size = int)) + geom_text(vjust = 2) +
  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed") +
  theme_bw() + xlab("") + ylab("") + ylim(c(-1, 12)) +
  theme(legend.position = "none")
dev.off()

m24 = module[[24]]
gk = myGK(m24)
m25 = module[[25]]
gk = myGK(m25)
data.frame(KEGG = gk$KEGG$Term[1:10], BP = gk$GO$BP$Term[1:10], MF = gk$GO$BP$Term[1:10])

# Master regulators with iRegulon
# !!!overlap with the WES variant genes!!!

# ---
load("bRNA/gene_expression_cnt.rdt") # DESeq2
load("bRNA/gene_expression_cnt_b6.rdt") # DESeq2
data <- gene_expression_cnt[rowMax(gene_expression_cnt) > 30, ] %>% as.data.frame
group <- factor(names(data), levels = c("BXSB", "BXSB_B6"))
colData <- data.frame(row.names = paste(names(data), 1:7, sep = "_"), condition = group)
countData <- as.integer(as.matrix(data))
countData <- t(apply(data, 1, as.integer))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resSig <- subset(res, padj < 0.05) %>% as.data.frame
plotMA(res, main = "DESeq2", ylim = c(-2,2))
summary(res)

# aligning to B6 and BXSB pseudo-genome does not make any difference
# ttests() on log2 TPM and DESeq() on counts returns very different gene lists
