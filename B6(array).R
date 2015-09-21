library(affy)
library(oligo)
library(dplyr)
library(genefilter)
library(preprocessCore)
library(Biobase)
library(ggplot2)
library(biomaRt)
library(qvalue)
library(quantro)
library(matrixStats)
library(Hmisc)
library(xlsx)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")
source("../../X/function.R")

# B6_YAA: MicroArray
cel.files <- read.celfiles(list.celfiles("bRNA/CHIP/CEL", full.name = T))
cel.rma <- rma(cel.files, target = "core") #
show(cel.rma)
b6 <- exprs(cel.rma)
colnames(b6) <- gsub("GC_Gene1-0ST_(.*)_430.*", "\\1", colnames(b6))

mogene <- read.delim("bRNA/CHIP/mogene2symbol.map", sep = "\t", stringsAsFactors = F)
b6 <- b6[rownames(b6) %in% mogene[, 1], ]
symbol <- mogene$Associated.Gene.Name[match(rownames(b6), mogene$Affy.MoGene.probeset)]
b6 <- apply(b6, 2, function(x) tapply(x, symbol, max))
b6 <- b6[rowMin(b6) > log2(120), ] %>% as.data.frame

sInf <- read.xlsx("bRNA/CHIP/Sample.xlsx", sheetName = "Sheet1", stringsAsFactors = F)
sInf <- sInf[match(names(b6), sInf$Sample.Name), c(1, 5, 6:7)]
sInf$Yaa = "Yaa"; sInf$Yaa[grep("wild type", sInf$Genotype)] = "B6"
sInf$Other = gsub(".*(Nfkb1|IL21R|TLR).*", "\\1", paste(sInf$Genotype, sInf$Other.Genotype))
sInf$Other[grep("(wild|Yaa)", sInf$Other)] = ""
sInf$Cell.Type <- gsub("^([BT]).*", "\\1", sInf$Cell.Type)
sInf$Group = paste(sInf$Yaa, sInf$Other, sep = "-")
sInf$Group = gsub("-$", "", sInf$Group)
rownames(sInf) = NULL

(summary = table(sInf$Cell.Type, sInf$Group))

pdf("bRNA/CHIP/sample.pdf")
barplot(summary, legend=T, beside=T, main="Sample Number by Geno- and Cell-type")
dev.off()

# marker genes: Tlr7, Nfkb1, Il21r, Tlr 
mutate(sInf, value = as.matrix(b6)["Tlr7", ]) %>% ggvis(~Group, ~value) %>% 
  layer_points(fill=~Cell.Type, shape=~Cell.Type) %>% layer_text(text:=~Sample.Name) %>% add_legend(c("fill", "shape"))

mutate(sInf, value = as.matrix(b6)["Nfkb1", ]) %>% ggvis(~Group, ~value) %>% 
  layer_points(fill=~Cell.Type) %>% layer_text(text:=~Sample.Name)

mutate(sInf, value = as.matrix(b6)["Il21r", ]) %>% ggvis(~Group, ~value) %>% 
  layer_points(fill=~Cell.Type) %>% layer_text(text:=~Sample.Name)

mutate(sInf, value = as.matrix(b6)["Tlr", ]) %>% ggvis(~Group, ~value) %>% 
  layer_points(fill=~Cell.Type) %>% layer_text(text:=~Sample.Name)

# cluster
rowVars(as.matrix(b6)) %>% summary

svd = svd(b6 - rowMeans(b6))
svd.pc <- diag(svd$d) %*% t(svd$v)

par(mfrow = c(1, 1))
plot(cumsum(svd$d) / sum(svd$d), ylim = c(0, 1), type = "b")

par(mfrow = c(2, 3))
for(x in 1:6) barplot(svd.pc[x, ], main = x)
for(x in 1:6) barplot(svd.pc[x, ], main = x, col = as.factor(sInf$Cell.Type))
legend("topright", levels(as.factor(sInf$Cell.Type)), fill=1:2, bg="white")
for(x in 1:6) barplot(svd.pc[x, ], main = x, col = as.factor(sInf$Group))
legend("topright", levels(as.factor(sInf$Group)), fill=1:5, bg="white")

km <- kmeans(t(b6), centers=10)
table(sInf$Group, km$cluster)
d <- dist(t(b6))
mds <- cmdscale(d)

mypar(1,2)
plot(mds[,1], mds[,2]) 
plot(mds[,1], mds[,2], col=km$cluster, pch=16)

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
