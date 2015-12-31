library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DNAcopy)
library(Gviz)

options(stringsAsFactors = F)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus/WES")

axTrack <- GenomeAxisTrack()
gTrackX <- IdeogramTrack(genome = "mm10", chromosome = "chrX")
gTrackY <- IdeogramTrack(genome = "mm10", chromosome = "chrY")

par <- list(chrX = c(169969759, 170931299), chrY = c(90745845, 91644698))
yaa <- list(chrX = c(161031299, 171031299), chrY = c(81644698, 91644698)) # 10Mbp

refGenes <- function(chr, from, to) UcscTrack(genome="mm10", chromosome=chr, track="refGene", from=from, to=to, 
  trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", name="Ensembl Genes", fill="red", 
  gene="name", symbol="name2", transcript="name", strand = "strand")

pdf("par.pdf", width = 12, height = 4)
txTrack <- refGenes("chrX", par$chrX[1], par$chrX[2])
plotTracks(list(gTrackX, axTrack, txTrack), from=par$chrX[1], to=par$chrX[2], transcriptAnnotation="symbol", stackHeight=0.2)
txTrack <- refGenes("chrY", par$chrY[1], par$chrY[2])
plotTracks(list(gTrackY, axTrack, txTrack), from=par$chrY[1], to=par$chrY[2], transcriptAnnotation="symbol", stackHeight=0.2)
dev.off()

pdf("yaa.pdf", width = 12, height = 4)

txTrack <- refGenes("chrX", yaa$chrX[1], yaa$chrX[2])
hlTrack <- HighlightTrack(trackList = txTrack, start = par$chrX[1], end = par$chrX[2], chromosome = "chrX", fill = "Skyblue")
plotTracks(list(gTrackX, axTrack, hlTrack), from=yaa$chrX[1], to=yaa$chrX[2], transcriptAnnotation="symbol", stackHeight=0.2)
symbol(txTrack) %>% unique

txTrack <- refGenes("chrY", yaa$chrY[1], yaa$chrY[2])
hlTrack <- HighlightTrack(trackList = txTrack, start = par$chrY[1], end = par$chrY[2], chromosome = "chrY", fill = "Skyblue")
plotTracks(list(gTrackY, axTrack, hlTrack), from=yaa$chrY[1], to=yaa$chrY[2], transcriptAnnotation="symbol", stackHeight=0.2)

dev.off()

pdf("chrY.pdf", width = 12, height = 4)
txTrack <- refGenes("chrY", 1, yaa$chrY[2])
hlTrack <- HighlightTrack(trackList = list(txTrack), start = par$chrY[1], end = par$chrY[2], chromosome = "chrY")
plotTracks(list(gTrackY, axTrack, hlTrack), from=1, to=yaa$chrY[2], transcriptAnnotation="symbol", stackHeight=0.2)
dev.off()

# variants and effects
sf_bxsb = "./BXSB_variants_filtered.vcf.gz"
scanVcfHeader(sf_bxsb)
vcf_bxsb = readVcf(sf_bxsb, genome = "mm10")
vcf_bxsb = vcf_bxsb[fixed(vcf_bxsb)$FILTER == "PASS"]
vcf_bxsb = vcf_bxsb[width(vcf_bxsb) == 1] # SNP only
vcf_bxsb = dropSeqlevels(vcf_bxsb, "chrM")
vcf_bxsb = renameSeqlevels(vcf_bxsb, gsub("chr", "", seqlevels(vcf_bxsb)))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdb <- renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))

all <- locateVariants(vcf_bxsb, txdb, AllVariants())

idx <- sapply(split(mcols(all)$QUERYID, mcols(all)$GENEID), unique)
variantByLocation <- sapply(names(idx), function(name) {
  d = all[mcols(all)$GENEID %in% name, c("QUERYID", "LOCATION")]
  table(mcols(d)$LOCATION[duplicated(d) == F])
})

rowSums(variantByLocation) 

sf = "/data/cgd/Sanger/REL-1505/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
(sanger_hdr = scanVcfHeader(sf))
ranges = GRanges(seqnames = 10, ranges = IRanges(start = 4455000, end = 4456000))
ranges = rowData(vcf_bxsb)
param = ScanVcfParam(geno = c("GT", "FI"), which = ranges)
vcf_sanger = readVcf(file = sf, genome = "mm10", param = param)
save(vcf_sanger, file = "./WES/sanger.rdt")
load("./WES/sanger.rdt")

sanger_geno <- geno(vcf_sanger)
sanger_geno$GT
sanger_geno$FI
		  
# !!! variants GWAS homolog, previous GWAS results (mouse/human)

# VARIANTS BY EXOME SEQ (BXSB)
var <- read.delim("./BXSB_SNP_filtered.vcf", skip = 23) # filter required
names(var)[1] <- "CHROM"

table(var$CHROM)

var$CHROM <- factor(var$CHROM, levels = c(1:19, "X", "Y"))

plot(table(var$CHROM), type = "h", lwd = 10, ylab = "Number")

chr <- read.delim("~/Dropbox/X/genomes/mouse.mm10.genome", sep = " ", stringsAsFactors = F, header = F)
chr$CHROM <- factor(gsub("chr", "", chr$V1), levels = c(1:19, "X", "Y"))

pdf("Yaa/genome.pdf", width = 18, height = 12)  # density per chromosome
ggplot() + # geom_density() +
  geom_segment(data = chr, aes(x = 1, y = 1e2, xend = V2, yend = 1e2), size = 5, color = "grey70") + 
  geom_bar(data = var, aes(x = POS), binwidth = 1e6, fill = "red") + 
  theme_bw() + xlab("") + ylab("") + facet_grid(CHROM ~ .) + 
  coord_cartesian(ylim=c(1, 999))
dev.off()

chrY <- filter(var, CHROM == "Y")  # vcf do not return YAA 

# EFFECTS PREDICTED BY ENSEMBL
vep <- read.delim("Yaa/BXSB_vep.txt", stringsAsFactors = F)
vep <- vep[apply(vep, 2, function(x) length(table(x)) != 1)]
vep$CHR <- gsub(":.*", "", vep$Location)
vep$POS <- gsub("^.*:(.*)-.*$", "\\1", vep$Location)  # CAVEAT INDELS

so <- read.delim("Yaa/so.txt", stringsAsFactors = F)
so.term <- gsub(" ", "", so$SO.term)  # variants by consequences
cons <- lapply(so.term, function(x) vep[grepl(x, vep$Consequence), ])
names(cons) <- so.term
cons <- cons[sapply(cons, function(x) nrow(x) > 0)]

(so.stat <- sort(sapply(cons, nrow), decreasing = T))

pdf("./Yaa/so.stat.pdf", width = 15, height = 10)
op <- par(mar = c(5, 20, 4, 2))
bar <- barplot(so.stat, xlim = c(0, 8e4), axes = F, border = NA, horiz = T, las = 1, space = 0.75)
abline(0, 1, lwd = 1, col = "black")
text(y = bar, x = so.stat + 2e3, labels = so.stat)
dev.off()

intron <- cons[["intron_variant"]]  # unique intron_variant, technical artifact
intron.stat <- tapply(intron$Consequence, intron$Location, function(x) unique(c(x))) %>% unlist %>% table
(intron.stat <- sort(intron.stat, decreasing = T))

# VarScan
CNA <- read.table("WES/BXSB_CNA.copynumber.called", header = T)
CNA.object <-CNA( genomdat = CNA[, 7], chrom = CNA[, 1], maploc = CNA[, 2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segment <- segment(CNA.smoothed, verbose = 1, min.width = 2)
segment <- segment(smoothed.CNA.object, undo.splits = "sdundo", undo.SD = 3,verbose = 1)

pdf("WES/segment.pdf", width = 10, height = 5)
plot(segment, plot.type="w")
plot(segment, plot.type="s")
dev.off()

plot(cov_chrX$V2, cov_chrX$V4, type = "l", ylim = c(0, 1e3))
plot(cov_chrX$V2, cov_chrX$V4, type = "l", xlim = yaa_loc$chrX)
