library(ape)
library(amap)
library(sva)
library(dplyr)
library(reshape)
library(ggvis)
library(ggplot2)
library(Biobase)
library(genefilter)
library(preprocessCore)
library(biomaRt)
library(quantro)
library(matrixStats)
library(xlsx)

rm(list = ls())
source("~/Dropbox/X/function.R")
setwd("~/Dropbox/GitHub/Lupus/B6(array)")

load("array.rdt") # Saved
for(obj in names(array)) assign(obj, array[[obj]])

# ---
library(affy); library(oligo)
cel.files <- read.celfiles(list.celfiles("CEL", full.name = T))
cel.rma <- rma(cel.files, target = "core"); show(cel.rma)
b6 <- exprs(cel.rma); colnames(b6) <- gsub("GC_Gene1-0ST_(.*)_430.*", "\\1", colnames(b6))

mogene <- read.delim("mogene2symbol.map", sep = "\t", stringsAsFactors = F)
b6 <- b6[rownames(b6) %in% mogene[, 1], ]
symbol <- mogene$Associated.Gene.Name[match(rownames(b6), mogene$Affy.MoGene.probeset)]
b6 <- apply(b6, 2, function(x) tapply(x, symbol, max))
b6 <- b6[rowMax(b6) > log2(120), ] %>% as.data.frame # Saved

sInf <- read.xlsx("Sample.xlsx", sheetName = "Sheet1", stringsAsFactors = F)
sInf <- sInf[match(names(b6), sInf$Sample.Name), c(1, 5, 6, 7)]

sInf$Yaa = "Yaa"; sInf$Yaa[grep("wild type", sInf$Genotype)] = "B6"
sInf$Other = gsub(".*(Nfkb1|IL21R|TLR).*", "\\1", paste(sInf$Genotype, sInf$Other.Genotype))
sInf$Other[grep("(wild|Yaa)", sInf$Other)] = ""
sInf$Cell.Type <- gsub("^([BT]).*", "\\1", sInf$Cell.Type)
sInf$Group = paste(sInf$Yaa, sInf$Other, sep = "-")
sInf$Group = gsub("-$", "", sInf$Group)
rownames(sInf) = sInf$Genotype = sInf$Other.Genotype = sInf$Yaa = sInf$Other = NULL

sInf_batch <- read.xlsx("Sample(DCR).xlsx", sheetName = "Sheet1", startRow = 2, stringsAsFactors = F)
sInf_batch <- sInf_batch[sInf_batch$Sample.Name %in% sInf$Sample.Name, ]
sInf$Sample.Prep.Sort.Date <- sInf_batch$Sample.prep.Sort.Date[match(sInf$Sample.Name, sInf_batch$Sample.Name)]
sInf$Affy.Date[sInf$Sample.Name %in% sInf_batch$Sample.Name] <- "20120409"
sInf$Sample.Prep.Sort.Date[29:48] = sInf$Affy.Date[29:48] = "20120530"

(summary = table(sInf$Cell.Type, sInf$Group))

pdf("sample1.pdf", width = 7, height = 5)
barplot(summary, legend=T, beside=T, main="Sample Number by Geno- and Cell-type")
dev.off()

table(sInf$Sample.Prep.Sort.Date, sInf$Affy.Date)

pdf("sample2.pdf", width = 7, height = 5)

ggplot(sInf, aes(x = Group, fill = Affy.Date)) + 
  geom_bar(width = 0.7) + facet_grid(~ Cell.Type) + theme_bw() +
  scale_fill_manual(values = c("dodgerblue3", "firebrick1")) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("") + ylab("Sample Number")

ggplot(sInf, aes(x = Group, fill = Sample.Prep.Sort.Date)) + 
  geom_bar(width = 0.7) + facet_grid(~ Cell.Type) + theme_bw() +
  scale_fill_manual(values = c("grey30", "firebrick1", "dodgerblue3", "chartreuse3")) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("") + ylab("Sample Number")

dev.off()

mutate(melt(t(b6)), group = sInf$Group[match(X1, sInf$Sample.Name)]) %>%
  ggplot(aes(X1,value,color=group)) + geom_boxplot() + xlab("") + theme(axis.text.x=element_text(angle=90))

.theme = function(data) 
  ggplot(data, aes(x = Group, y = value, color = Group, label = Sample.Name)) +
  geom_point(position = position_jitter(w = 0.2, h = 0.0)) + geom_text(size = 1, color = "grey30") +
  facet_grid(Cell.Type ~ Sample.Prep.Sort.Date) + xlab("") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) 

pdf("genotype.pdf", width = 7, height = 5)

mutate(sInf, value = as.matrix(b6)["Tlr7", ]) %>% .theme()
mutate(sInf, value = as.matrix(b6)["Nfkb1", ]) %>% .theme()
mutate(sInf, value = as.matrix(b6)["Il21r", ]) %>% .theme()

dev.off()

# Labeling and genotyping issues
b6 = b6[sInf$Sample.Prep.Sort.Date != "20120530"]
sInf = sInf[sInf$Sample.Prep.Sort.Date != "20120530", ]

b6_t = b6[sInf$Cell.Type == "T"]
sInf_t = sInf[sInf$Cell.Type == "T", ]

# Cluster

svd = svd(b6 - rowMeans(b6))
svd = svd(b6_t - rowMeans(b6_t))

svd.pc <- diag(svd$d) %*% t(svd$v)

pdf("../pdf/batch.pdf", width = 8, height = 5) # all samples

par(mfrow = c(1, 2))
plot(svd$d^2 / sum(svd$d^2), ylim = c(0, 1), type = "b", ylab = "Percentage variation")
boxplot(split(svd$v[, 1], as.factor(sInf$Cell.Type)))

plot(t(cor(as.numeric(as.factor(sInf$Sample.Prep.Sort.Date)), svd$v)), ylab = "Correlation")

cols = as.numeric(as.factor(sInf$Sample.Prep.Sort.Date))
plot(svd$v[, 2], svd$v[, 3], col = cols, pch = 16, xlab = "PC 2", ylab = "PC 3")
legend("bottomright", levels(as.factor(sInf$Sample.Prep.Sort.Date)), fill=1:3, cex = 0.7, bg="white")

dev.off()

# B/T so different, would this swamp genotype signal?

pdf("../pdf/batch2.pdf", width = 8, height = 5) # t cell only

par(mfrow = c(1, 2))
plot(svd$d^2 / sum(svd$d^2), ylim = c(0, 1), type = "b", ylab = "Percentage variation")

plot(t(cor(as.numeric(as.factor(sInf_t$Sample.Prep.Sort.Date)), svd$v)), ylab = "Correlation")
text(t(cor(as.numeric(as.factor(sInf_t$Sample.Prep.Sort.Date)), svd$v)), labels = 1:14, pos = 1)

par(mfrow = c(1, 1))
cols = as.numeric(as.factor(sInf_t$Sample.Prep.Sort.Date))
plot(svd$v[, 1], svd$v[, 2], col = cols, pch = 16, xlab = "PC 1", ylab = "PC 2")
legend("bottomright", levels(as.factor(sInf_t$Sample.Prep.Sort.Date)), fill=1:3, cex = 0.7, bg="white")

dev.off()

par(mfrow = c(2, 3))
for(x in 1:6) barplot(svd.pc[x, ], main = x, col = as.factor(sInf$Cell.Type))
legend("topright", levels(as.factor(sInf$Cell.Type)), fill=1:2, bg="white")

for(x in 1:6) barplot(svd.pc[x, ], main = x, col = as.factor(sInf$Group))
legend("topright", levels(as.factor(sInf$Group)), fill=1:5, bg="white")

# Model
(cell = factor(sInf$Cell.Type, levels = c("B", "T")))
(genoAll = factor(sInf$Group, levels = c("B6", "Yaa", "Yaa-IL21R", "Yaa-Nfkb1")))
(genoYaa = factor(gsub("(B6|Yaa).*", "\\1", sInf$Group), levels = c("B6", "Yaa")))

genoIL21R = factor(rep("WT", 28), levels = c("WT", "IL21R")); genoIL21R[grep("IL21R", sInf$Group)] = "IL21R"
genoNfkb1 = factor(rep("WT", 28), levels = c("WT", "Nfkb1")); genoNfkb1[grep("Nfkb1", sInf$Group)] = "Nfkb1"

(batchSort = factor(sInf$Sample.Prep.Sort.Date))
(batchAffy = factor(sInf$Affy.Date))

# SVA 
(model = model.matrix(~ cell + cell * genoYaa + cell * genoIL21R + cell * genoNfkb1))
sva_fit = sva(as.matrix(b6), model)

head(sva_fit$sv)
cor(sva_fit$sv[, 1], svd$v)
cor(sva_fit$sv[, 2], svd$v)
cor(sva_fit$sv[, 3], svd$v)

(X = model.matrix(~ cell + cell * genoYaa + cell * genoIL21R + cell * genoNfkb1 + sva_fit$sv))

fit0 <- apply(b6, 1, function (x) summary(lm(x ~ X - 1)))

# Use residules to cluster samples

sv = sapply(fit0, function(x) x$coefficients[c("Xsva_fit$sv1", "Xsva_fit$sv2", "Xsva_fit$sv3"), "Estimate"])
b6_res = as.matrix(b6) - t(sva_fit$sv %*% sv) 

mycol = as.numeric(as.factor(sInf$Sample.Prep.Sort.Date))
hc1 <- hcluster(t(b6), method = "pearson", link = "average") %>% as.phylo
hc2 <- hcluster(t(b6_res), method = "pearson", link = "average") %>% as.phylo
hc1$tip.label = hc2$tip.label = paste(sInf$Cell.Type, sInf$Group, sep = "-")

pdf("../pdf/hc1.pdf", width = 12, height = 5)

par(mar = c(5, 0, 4, 2), mfrow = c(1, 1))
plot(hc1, cex=0.7, tip.color = mycol)
legend("top", levels(as.factor(sInf$Sample.Prep.Sort.Date)), cex = 0.7, fill=1:3, bg="white")

plot(hc2, cex=0.7, tip.color = mycol)
legend("top", levels(as.factor(sInf$Sample.Prep.Sort.Date)), cex = 0.7, fill=1:3, bg="white")

dev.off()

# decompose causal effects

(X = model.matrix(~ cell + cell * genoYaa + cell * genoIL21R + cell * genoNfkb1))

fit1 <- apply(b6_res, 1, function (x) summary(lm(x ~ X - 1)))

fit.e = t(sapply(fit1, function(x) x$coefficients[-1, "Estimate"]))
fit.p = t(sapply(fit1, function(x) x$coefficients[-1, "Pr(>|t|)"]))

sig = apply(fit.p, 2, function(x) sum(x < 0.05))
names(sig) = gsub("^X", "", names(sig))

pdf("../pdf/sig_stat.pdf", width = 8, height = 5)

op <- par(mar = c(10, 4, 4, 2))
bar <- barplot(sig, ylim = c(0, 1e4), axes = F, border = NA, las = 2, col = "blue")
abline(0, 1, lwd = 1, col = "black")
text(x = bar, y = sig + 1e1, labels = sig)

dev.off()

# genes and KEGG

(sig_b = names(which(fit.e[, 1] < 0 & fit.p[, 1] < 0.05)))
(sig_t = names(which(fit.e[, 1] > 0 & fit.p[, 1] < 0.05)))

(sig_t_Nfkb1 = names(which(fit.p[, 7] < 0.05)))

gk = gk_b = mmGK(sig_b)
gk = gk_t = mmGK(sig_t)
gk = gk_t_Nfkb1 = mmGK(sig_t_Nfkb1)
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

gk_t = gk

# batch correction significantly increase adjusted R-squared??

# adjusted R2 keep falling

lm(x ~ cell + genoYaa + genoIL21R + genoNfkb1) %>% summary
lm(x ~ cell + genoYaa + genoTLR + genoIL21R + genoNfkb1 + cell * genoTLR)
lm(x ~ cell + genoYaa + genoTLR + genoIL21R + genoNfkb1 + cell * genoTLR + cell * genoIL21R + cell * genoNfkb1)
lm(x ~ cell + genoYaa + genoTLR + cell * genoTLR) # Now, IL21R and Nfkb1 taken as Yaa

geneAll = rownames(fit.p)
myList = sapply(colnames(fit.p), function(x) geneAll[fit.p[, x] < 0.05])

myList = sapply(colnames(fit.p), function(x) {
  y = cbind(Pval = fit.p[, x], Size = fit.e[, x])
  y[y[, "Pval"] < 0.05, ]
})

names(myList) <- gsub(":", "-", names(myList))

lapply(myList, nrow)

lapply(names(myList), function(x) write.xlsx(myList[[x]], file = "../myfile.xlsx", sheetName = x, append = T))

b6_avg = sapply(group_all, function(x) rowMeans(b6_res[, sInf$Group == x]))

group_all = c( "B6", "Yaa", "Yaa-Nfkb1", "Yaa-IL21R")
write.xlsx(b6_avg, file = "../myfile.xlsx", sheetName = "Expression", append = T)

