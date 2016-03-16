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

rm(list = ls())
source("~/Dropbox/GitHub/X/function.R")
setwd("~/Dropbox/GitHub/Lupus/BXSB(rnaseq)")

load("gene_expression_tpm.rdt"); tpm = gene_expression_tpm # aligned to BXSB
load("gene_expression_tpm_b6.rdt"); tpm_b6 = gene_expression_tpm # aligned to B6

rm(gene_expression_tpm)
select = intersect(rownames(tpm), rownames(tpm_b6))

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mart1 = getBM(c("ensembl_gene_id", "external_gene_name"), "ensembl_gene_id", select, mart)

select = select[select %in% mart1$ensembl_gene_id]
tpm = tpm[select, ]
tpm_b6 = tpm_b6[select, ]

rownames(tpm) = rownames(tpm_b6) = mart1$external_gene_name

# predicted genes show up: KB's pipeline

x = rowMeans(log2(tpm + 1))
y = rowMeans(log2(tpm_b6 + 1))
z = cbind(sb = x, b6 = y)

z = as.data.frame(z[rowMax(z) > 3, ])

# MD plot

plot(rowSums(z), z$sb-z$b6) # X: sum; Y: log fold

y1 = (z$sb - z$b6) / rowSums(z) # log fold normalized by sum
plot( rowSums(z), y1 )

table(abs(y1) == 1)
table(abs(y1) > 0.9)

points(rowSums(z)[abs(y1) > 0.9], y1[abs(y1) > 0.9], col = "red", pch = 8)

plot(rowSums(z), z$sb-z$b6) # X: sum; Y: log fold
points(rowSums(z)[abs(y1) > 0.9], (z$sb-z$b6)[abs(y1) > 0.9], col = "red", pch = 8)

(z1 = z[abs(y1) > 0.9, ])

(genes = rownames(z1)[! grepl("^Gm", rownames(z1))])

myGk = mmGK(genes)
