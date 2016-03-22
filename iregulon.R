library(igraph)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")
options(stringsAsFactors = F)

# PARSE THE IREGULON OUTPUT

file <- "./iregulon/gene1.tsv"
univ <- read.table("./iregulon/universe.txt")$V1

# about univ: all genes with maximal expression level at least 10 in tpm, refer to emase.R

ireg <- read.delim(file, comment.char = ";")
factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))

factor <- lapply(factor, function(x) x[x %in% univ])
idx <- sapply(factor, length) > 0

factor <- factor[idx]; target <- target[idx]

edges <- lapply(1:length(factor), function(x) expand.grid(factor[[x]], target[[x]], stringsAsFactors = F))
edges <- do.call(rbind, edges)
edges <- edges[! duplicated(edges), ]

# VISUALIZATION: IGRAPH
igraph.dt <- graph.data.frame(edges)

igraph.dt$layout <- layout.sphere
igraph.dt$layout <- layout.circle
igraph.dt$layout <- layout.fruchterman.reingold

V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
V(igraph.dt)$color[V(igraph.dt)$name %in% unlist(factor)] <- "gold"
V(igraph.dt)$size = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 10 + 2
V(igraph.dt)$label.cex = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 2 + 0.5
V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"

pdf("./iregulon/gene1.pdf", width = 9, height = 9)

plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)

dev.off()
