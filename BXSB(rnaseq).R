library(dplyr)
library(Biobase)
library(genefilter)
library(preprocessCore)
library(ggplot2)
library(biomaRt)
library(quantro)
library(qvalue)
library(matrixStats)
library(Hmisc)
library(ape)
library(amap)

rm(list = ls())
source("~/Dropbox/GitHub/X/function.R")
setwd("~/Dropbox/GitHub/Lupus/BXSB(rnaseq)")

load("./gene_expression_tpm.rdt")
tpm = gene_expression_tpm # aligned to BXSB
rm(gene_expression_tpm)

tpm = tpm[rowMax(tpm) > 8, ]

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mart1 = getBM(c("ensembl_gene_id", "external_gene_name"), "ensembl_gene_id", rownames(tpm), mart)

tpm = tpm[mart1$ensembl_gene_id, ]
rownames(tpm) = mart1$external_gene_name

data <- tpm[rowMax(tpm) > 10, ] %>% as.data.frame
(group <- factor(names(data), levels = c("BXSB", "BXSB_B6")))

outlier <- data[rowMax(as.matrix(data)) > 1e4, ]

apply(log2(outlier + 1), 1, function(x) summary(lm(x ~ group))$coefficients["groupBXSB_B6", "Pr(>|t|)"])

# Igkc and Ighg2c are two immunity genes that showed dramatic change between the two strains

data <- data[! rowMax(as.matrix(data)) > 1e4, ]
data <- sweep(data, 2, colSums(data), "/") * 1e6

pdf("bRNA/data.pdf", width = 8, height = 6)
par(mar = c(10, 4, 4, 2))
matboxplot(log2(data + 1), groupFactor = group)
dev.off()

# Clustering

hc1 <- hcluster(t(data), method = "pearson", link = "average")
plot(as.phylo(hc1), type = "fan")

# Pair-wise ttest

ttest <- rowttests(as.matrix(log2(data + 1)), group) # BXSB/B6 -> BXSB

data$BXSB_B6_avg <- rowMeans(data[group == "BXSB_B6"])
data$BXSB_avg <- rowMeans(data[group == "BXSB"])

data <- cbind(data, ttest)
data$q.value <- qvalue(data$p.value)$qvalues

data$sel <- "N"
data$sel[abs(data$dm) > 0.3 & data$q.value < 0.05] <- "S"
# !!! transcriptome of BXSB_B6 and BXSB so different!

pdf("bRNA/vocano.pdf", width = 6, height = 4)

ggplot(data, aes(x = dm, y = -log10(p.value))) +
  geom_point(aes(color = as.factor(sel)), size = 1) +
  scale_color_manual(values = c("grey30", "firebrick1")) +
  theme_bw() + xlab("Effect") + ylab("-log10(pvalue)") +
  theme(legend.position = "none")

dev.off()

genes = rownames(data)[data$sel == "S"]

myGk <- mmGK(genes)

myGk$KEGG[1:10, ]
sapply(myGk$GO, function(x) x$Term[1:20])
data.frame(KEGG = myGk$KEGG$Term[1:20], BP = myGk$GO$BP$Term[1:20], MF = myGk$GO$BP$Term[1:20])

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Median norm.

objectMedians <- colMedians(as.matrix(data))
normalize.median <- sweep(data, 2, objectMedians, FUN = "-")
colMedians(as.matrix(normalize.median))

# Quantile norm.

normalize.quantile <- normalize.quantiles(as.matrix(data), copy = TRUE)
dimnames(normalize.quantile) <- dimnames(data)

par(mfrow = c(1, 2))
qqplot(data[, 1], data[, 7], xlab = "BXSB_B6", ylab = "BXSB")
abline(0, 1)
qqplot(normalize.median[, 1], normalize.median[, 7], xlab = "BXSB_B6", ylab = "BXSB")
abline(0, 1)
qqplot(normalize.quantile[, 1], normalize.quantile[, 7], xlab = "BXSB_B6", ylab = "BXSB")
abline(0, 1)

