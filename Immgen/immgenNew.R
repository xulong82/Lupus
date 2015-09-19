library(dplyr)
library(ggplot2)
library(xlsx)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

file <- "Immgen/RequestedImmGenData2015-08-11_13-20-24.csv.gz"
immgen <- read.csv(file, stringsAsFactors = F, header = T, check.names = F) 
immgen <- apply(immgen[4:227], 2, function(x) tapply(x, immgen$GeneSymbol, max))
immgen <- log2(immgen + 1) %>% as.data.frame

(type <- gsub("_.*", "", names(immgen)))
(tissue <- gsub(".*_", "", names(immgen)))

module <- read.xlsx("Immgen/gene_assignment.xls", sheetIndex = 1, stringsAsFactors = F)
table(module$Coarse.module)
module <- lapply(1:81, function(x) module$Gene[module$Coarse.module == x])
names(module) <- paste0("M", 1:81)

immgenList = list(immgen = immgen, module = module)
save(immgenList, file = "./Immgen/immgenList.rdt")

# ---

rnaseq <- myTpm[geneId, ]
rnaseq <- log2(rnaseq + 1)
rnaseq <- rnaseq[apply(rnaseq, 1, max) > 5, ]
rnaseq <- data.frame("NN"=rowMeans(rnaseq[1:2]), "NP"=rowMeans(rnaseq[3:4]), "PP"=rowMeans(rnaseq[5:6]))

idx <- intersect(rownames(rnaseq), rownames(immgen))  
data1 <- cbind(rnaseq[idx, ], immgen[idx, ])
data1 <- as.data.frame(apply(data1, 2, scale))  
rownames(data1) <- idx

pearson <- cor(data1, method = "pearson")  
heatmap(pearson[1:3, ], cexRow = 1, cexCol = .2)

pearson_NN20 <- data.frame(Sample = names(sort(pearson[, "NN"], decreasing = T)[1:20]))
pearson_NP20 <- data.frame(Sample = names(sort(pearson[, "NP"], decreasing = T)[1:20]))
pearson_PP20 <- data.frame(Sample = names(sort(pearson[, "PP"], decreasing = T)[1:20]))

data2 <- data1[, union(pearson_NN20$Sample, union(pearson_NP20$Sample, pearson_PP20$Sample))]
data2 <- cbind(data2[c("NN","NP","PP")], data2[grep("_", colnames(data2))])
colnames(data2)[1:3] <- c("N", "ACT", "IL21")
heatmap(cor(data2))

pdf("Results/immgenLevel.pdf", width = 5, height = 8)
levelplot(cor(data2)[1:3, -c(1:3)], main = "", xlab   = "", ylab = "")
dev.off()

#--- IMMGEN AND RNASEQ PCA ---
sample1 <- colnames(data2)[grep("T_4Nve_", colnames(data2))]
sample2 <- colnames(data2)[grep("NKT_4_", colnames(data2))]
sample3 <- colnames(data2)[grep("T_4Mem", colnames(data2))]

data3 <- data2[c("N", "ACT", "IL21", sample1, sample2, sample3)]
pca <- prcomp(t(data3), scale = TRUE)  # PCA

pca.x <- pca$x[, 2:3]
type <- c("NN", "NP", "PP", rep("T_Nve", 4), rep("T_NK", 4), rep("T_Mem", 4))
immgen.pca <- data.frame(pca.x, type)
save(immgen.pca, file = "./data/immgenPCA.rdt")

par(mfrow = c(2, 3), mar = c(5, 4, 4, 2), las = 3)
for (i in 1:6) {
  name <- paste("PCA", i, sep = " ")
  barplot(pca$x[, i],  main = name) 
  abline(0, 0, lwd = 1, col = "black")
}

# --- GRAPH ---
load("./data/immgenPCA.rdt")
immgen.pca$lab <- rep(1, nrow(immgen.pca))
immgen.pca$lab[1:3] <- 2
mycolors <- c("darkorchid4", "darkgreen", "darkblue", "dodgerblue", "tan4", "brown2")
immgen.pca$PC2 <- -immgen.pca$PC2
immgen.pca$PC3 <- -immgen.pca$PC3
pdf("Results/immgenPCA.pdf", fonts = "Helvetica", height = 7, width = 8)
ggplot(immgen.pca, aes(x = PC3, y = PC2)) + 
  geom_point(aes(color = type, shape = as.factor(lab), size = lab), scale_size = 10) +
  xlim(-26, 26) + ylim(-26, 26) + theme_bw() +
  scale_size_continuous(range = c(5,10), guide = FALSE) +
  scale_shape_discrete(guide = FALSE) +
  scale_color_manual(values = mycolors,
                     breaks = c("NN", "NP", "PP", "T_Mem", "T_NK", "T_Nve"), 
                     labels = c("N", "ACT", "ACT IL21", "T Memory", "NKT", "Naive T")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 1, color = "black"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(vjust = -0.5),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
