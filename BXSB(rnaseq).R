library(dplyr)
library(Biobase)
library(genefilter)
library(preprocessCore)
library(ggplot2)
library(biomaRt)
library(qvalue)
library(quantro)
library(matrixStats)
library(Hmisc)

rm(list = ls())
source("~/Dropbox/X/function.R")
setwd("~/Dropbox/GitHub/Lupus/BXSB(rnaseq)")

load("gene_expression_tpm.rdt"); tpm = gene_expression_tpm
load("gene_expression_tpm_b6.rdt"); tpm_b6 = gene_expression_tpm

select = intersect(rownames(tpm), rownames(tpm_b6))

# predicted genes show up: KB's pipeline
x = log2(tpm[select, 4] + 1)
y = log2(tpm_b6[select, 4] + 1)
plot( (x+y)/2, x-y)

data <- gene_expression_tpm[rowMax(gene_expression_tpm) > 10, ] %>% as.data.frame
(group <- factor(names(data), levels = c("BXSB", "BXSB_B6")))

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
