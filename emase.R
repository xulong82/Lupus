library(dplyr)
library(xlsx)
library(biomaRt)
library(Biobase)
library(genefilter)
library(VennDiagram)
library(qvalue)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")
source("~/Dropbox/GitHub/X/function.R")
options(stringsAsFactors = F)

# Rna-seq on bxsb and bxsb/b6y

name <- list.files(path = "./emase/", pattern = "BXSB_.*.tpm")
rna.emase <- lapply(name, function(x) {
  cat(x, "\n"); filepath = file.path("./emase/", x); read.delim(filepath)
}); names(rna.emase) <- gsub("_GES.*", "", name)

rna.emase.B6 <- sapply(rna.emase, function(x) x$B) %>% as.data.frame
rna.emase.SB <- sapply(rna.emase, function(x) x$S) %>% as.data.frame
rownames(rna.emase.B6) <- rownames(rna.emase.SB) <- rna.emase[[1]]$locus

# Clean the data

data = cbind(rna.emase.B6, rna.emase.SB)
data = data[rowMax(as.matrix(data)) > 10, ]

mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(c("ensembl_gene_id", "external_gene_name"), "ensembl_gene_id", rownames(data), mart)

data = data[biomart$ensembl_gene_id, ]
rownames(data) = biomart$external_gene_name

yaa = factor(rep(c(rep("yaa0", 3), rep("yaa1", 4)), 2), levels = c("yaa0", "yaa1"))
str = factor(c(rep("B6", 7), rep("SB", 7)), levels = c("B6", "SB"))

myplot = function(data, geneId) {
  x = data.frame(value = as.matrix(data)[geneId, ], yaa = yaa, str = str)
  ggplot(x, aes(x = yaa, y = value, fill = str)) +
  geom_boxplot() + theme_bw() + xlab("") + ggtitle(geneId) +
  scale_fill_manual(values = c("white", "firebrick1"))
}

pdf("./pdf/example.pdf", width = 3, height = 4)

myplot(data, "Tlr7") # B6 genes with Yaa effects
myplot(data, "Tlr2") # SB genes with Yaa effects
myplot(data, "Tlr9") # SB genes without Yaa effects

dev.off()

# define the SB genes

b6 = rowMeans(data[str == "B6"])
sb = rowMeans(data[str == "SB"])

pdf("./pdf/gene.pdf", width = 6, height = 4)

par(mfrow = c(1, 2))
plot(b6 + sb, b6 - sb, xlim = c(0, 1e4))
abline(0, 1, col = "red"); abline(0, -1, col = "red")
plot(b6 + sb, b6 - sb, xlim = c(0, 1e3), ylim = c(-1e3, 1e3))
abline(0, 1, col = "red"); abline(0, -1, col = "red")

dev.off()

# mean(data_sb) + mean(data_b6) = 1,

data2 = data / rowMeans(data) / 2

ttest = rowttests(as.matrix(data2), str)
ttest$q.value = p.adjust(ttest$p.value, method = "fdr")

# dm as 0.9 tells more than 95% of the reads are strain-specific

gene.sb = rownames(ttest)[ttest$q.value < 0.05 & (ttest$dm < -0.9)]
gene.b6 = rownames(ttest)[ttest$q.value < 0.05 & (ttest$dm > 0.9)]

vennList <- list(SB = gene.sb, B6 = gene.b6)
venn.diagram(vennList, imagetype = "png", file = "./pdf/venn.png", width = 2000, height = 2000)

myGk = mmGK(gene.sb)
myGk = mmGK(gene.b6)
data.frame(KEGG = myGk$KEGG$Term[1:20], BP = myGk$GO$BP$Term[1:20])

# define the YAA genes

data3 = data[1:7] + data[8:14]
ttest = rowttests(as.matrix(log2(data3 + 1)), yaa[1:7])
ttest$q.value = p.adjust(ttest$p.value, method = "fdr")

# dm is log2 fold change, 0.1 equals a 7% change

gene.yaa = rownames(ttest)[ttest$q.value < 0.05 & (abs(ttest$dm) > 0.1)]

# venn

vennList <- list(SB = gene.sb, B6 = gene.b6, YAA = gene.yaa)
venn.diagram(vennList, imagetype = "png", file = "./pdf/venn.png", cat.pos = 0, width = 2000, height = 2000)

myGk = mmGK(intersect(gene.sb, gene.yaa))
data.frame(KEGG = myGk$KEGG$Term[1:20], BP = myGk$GO$BP$Term[1:20], MF = myGk$GO$BP$Term[1:20])

write(intersect(gene.sb, gene.yaa), file = "./iregulon/gene1.txt")

# B6 and SB genes produce different proteins
# yaa genes orchestrate protein amounts

# wes on bxsb strain

wes.emase <- read.delim("./emase/emase.isoforms.effective_read_counts")
wes.emase <- wes.emase[rowMax(as.matrix(wes.emase[c("B", "S")])) > 50, ]
wes.emase <- data.frame(row.names = wes.emase$locus, B6 = wes.emase$B, SB = wes.emase$S)

mart1 <- getBM(c("ensembl_gene_id", "external_gene_name"), "ensembl_gene_id", rownames(wes.emase), mart)
mart1 <- mart1[! duplicated(mart1$external_gene_name), ]

wes.emase <- wes.emase[mart1$ensembl_gene_id, ]
rownames(wes.emase) <- mart1$external_gene_name

wes.emase$diff <- (wes.emase$S - wes.emase$B) / rowMeans(wes.emase) / 2

pdf("./pdf/hist.pdf", width = 6, height = 4)

hist(wes.emase$diff)

dev.off()

wes.sb.gene = rownames(wes.emase)[wes.emase$diff > 0.9]
wes.b6.gene = rownames(wes.emase)[wes.emase$diff < -0.9]

vennList <- list(RNA_SB = gene.sb, RNA_B6 = gene.b6, WES_SB = wes.sb.gene, WES_B6 = wes.b6.gene)
venn.diagram(vennList, imagetype = "png", file = "./pdf/venn.png", cat.pos = 0, width = 2000, height = 2000)

# because the predictors are dependent, B6/yaa0 and B6/yaa1, SB/yaa0 and SB/yaa1, linear regression is not appropriate

# fit <- apply(data2, 1, function (x) lm(x ~ str + yaa + str * yaa))
# fit <- apply(log2(data + 1), 1, function (x) lm(x ~ str + yaa + str * yaa))

# summary(fit[["Tlr7"]])
#
# fit.fs <- sapply(fit, function (x) summary(x)$fstatistic)
# fit.pv <- apply(fit.fs, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F))
# fit1 <- fit[fit.pv < 0.05]
#
# fit1.e = t(sapply(fit1, function(x) summary(x)$coefficients[-1, "Estimate"]))
# fit1.p = t(sapply(fit1, function(x) summary(x)$coefficients[-1, "Pr(>|t|)"]))
#
# logit = ( fit1.p < 0.05 )
# logit = ( fit1.p < 0.05 & abs(fit1.e) > 0.3 )
# logit = logit[rowSums(logit) != 0, ]
#
# profile <- apply(logit, 1, function (x) paste(x, collapse = "-"))
# (profile.table <- sort(table(profile)))
# profile.name <- names(profile.table)
# profile.gene = sapply(profile.name, function(x) names(profile)[profile == x])
#
# tile.dt = do.call(rbind, lapply(profile.name, function(x) as.logical(unlist(strsplit(x, "-")))))
# tile.dt <- data.frame(v = c(tile.dt), p = rep(profile.name, 3), g = rep(colnames(logit), each = nrow(tile.dt)))
# tile.dt$p = factor(tile.dt$p, levels = profile.name)
# tile.dt$g = factor(tile.dt$g, levels = colnames(logit))
#
# ggplot(tile.dt, aes(x = g , y = p , fill = v)) + geom_tile(colour = "white") +
#   theme_bw() + xlab("") + ylab("") + coord_flip() +
#   # scale_x_discrete(labels = c("2m:APP", "4m:APP", "5m:APP")) +
#   scale_y_discrete(labels = profile.table) +
#   scale_fill_manual(values = c("grey80", "firebrick1"))
