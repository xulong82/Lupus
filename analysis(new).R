
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

# @@@@@@@@@@@@@@@@@@@@@@@@

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
