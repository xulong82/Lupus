rm(list = ls())

setwd("~/Dropbox/GitHub/Lupus")

fname <- list.files(path = "/data/xwang/Lupus/RSEM_BXSB", pattern = "*.genes.results")
fname <- list.files(path = "/data/xwang/Lupus/RSEM_B6", pattern = "*.genes.results")

gene_expression <- lapply(fname, function(x) { cat(x, "\n")
# filepath = file.path("/data/xwang/Lupus/RSEM_B6", x)
  filepath = file.path("/data/xwang/Lupus/RSEM_BXSB", x)
  read.delim(filepath, stringsAsFactors = F)[, -2]
}); names(gene_expression) <- gsub("_Y.*", "", fname)

gene_expression_tpm <- sapply(gene_expression, function(x) x$TPM)
rownames(gene_expression_tpm) <- gene_expression[[1]]$gene_id

gene_expression_cnt <- sapply(gene_expression, function(x) x$expected_count)
rownames(gene_expression_cnt) <- gene_expression[[1]]$gene_id

setwd("~/Dropbox/GitHub/Lupus/BXSB(rnaseq)")

save(gene_expression_tpm, file = "gene_expression_tpm.rdt")
save(gene_expression_cnt, file = "gene_expression_cnt.rdt")
save(gene_expression_tpm, file = "gene_expression_tpm_b6.rdt")
save(gene_expression_cnt, file = "gene_expression_cnt_b6.rdt")

list.files()
sample  <- read.csv("trim.sh.o622832", stringsAsFactors = F, header = F)
sample <- gsub("_GES.*", "", sample$V1[2:8])

trim <- read.csv("trim.sh.e622832", stringsAsFactors = F, header = F)
trim <- trim$V1[grep("^Input", trim$V1)]
trim <- gsub(".*\\(", "", gsub("%.*","", trim)) %>% as.numeric

files <- list.files("log_rsem_b6", pattern = ".sh.e")
align_b6 <- lapply(files, function(x) read.csv(paste0("log_rsem_b6/", x), stringsAsFactors = F, header = F))
align_b6 <- lapply(align_b6, function(x) x$V1[grep("at least one", x$V1)])
align_b6 <- sapply(align_b6, function(x) gsub(".*\\((.*)%.*", "\\1", x)) %>% as.numeric

files <- list.files("log_rsem_pseudo_bxsb", pattern = ".sh.e")
align_bxsb <- lapply(files, function(x) read.csv(paste0("log_rsem_pseudo_bxsb/", x), stringsAsFactors = F, header = F))
align_bxsb <- lapply(align_bxsb, function(x) x$V1[grep("at least one", x$V1)])
align_bxsb <- sapply(align_bxsb, function(x) gsub(".*\\((.*)%.*", "\\1", x)) %>% as.numeric

df <- data.frame(Sample = sample, Trim = trim, B6 = align_b6, BXSB = align_bxsb)
df <- df[c(4, 2, 5, 7, 3, 6, 1), ]

pdf("pipeline.pdf", width = 8, height = 5)
plot(df$Trim, type = "b", col = "black", xaxt = "n", xlab = "", ylab = "%", ylim = c(70, 100))
lines(df$B6, type = "b", col = "blue")
lines(df$BXSB, type = "b", col = "red")
legend("topright", c("Trim","Align(B6)","Align(BXSB)"), fill = c("black", "blue", "red"), horiz = TRUE)
axis(1, at = 1:7, labels = gsub("BXSB_", "", df$Sample))
dev.off()
