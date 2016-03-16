library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ensemblVEP)

options(stringsAsFactors = F)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
seqlevels(txdb) <- gsub("chr", "", seqlevels(txdb))

# compare the sanger and jax/wes variants of bxsb

sf = "./wes/BXSB_variants_filtered.vcf.gz"
sf = "./Sanger/Sanger.UNC.Combined.SNPs.mm10.BXSB.vcf.gz"

scanVcfHeader(sf)
vcf = readVcf(sf, genome = "mm10")
seqlevels(vcf)

vcf = vcf[fixed(vcf)$FILTER == "PASS"]
vcf = vcf[width(vcf) == 1] # SNP only
vcf = dropSeqlevels(vcf, "chrM")
vcf = renameSeqlevels(vcf, gsub("chr", "", seqlevels(vcf)))

all <- locateVariants(vcf, txdb, AllVariants())

all.wes= all
all.sanger = all

table(mcols(all.wes)$LOCATION)
table(mcols(all.sanger)$LOCATION)

# coding::overlaps

all.wes = all.wes[width(all.wes) == 1] # SNP only

coding.wes = all.wes[mcols(all.wes)$LOCATION == "coding"]
coding.sanger = all.sanger[mcols(all.sanger)$LOCATION == "coding"]

length(unique(mcols(coding.wes)$QUERYID))
length(unique(mcols(coding.sanger)$QUERYID))

xx = intersect(coding.wes, coding.sanger)

length(xx)
length(unique(xx))

length(coding.wes)
length(unique(coding.wes))

length(coding.sanger)
length(unique(coding.sanger))

length(xx) / length(unique(coding.wes))
length(xx) / length(unique(coding.sanger))

# Note: Sanger/wes have large overlap

# Next: Sanger's SB and WES's BXSB

# Ensembl VEP

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Backup

var <- read.delim("./BXSB_SNP_filtered.vcf", skip = 23) # filter required
names(var)[1] <- "CHROM"

table(var$CHROM)
var$CHROM <- gsub("chr", "", var$CHROM)
var$CHROM <- factor(var$CHROM, levels = c(1:19, "X", "Y"))

plot(table(var$CHROM), type = "h", lwd = 10, ylab = "Number")

chr <- read.delim("~/Dropbox/GitHub/X/genomes/mouse.mm10.genome", sep = " ", header = F)
chr$CHROM <- factor(gsub("chr", "", chr$V1), levels = c(1:19, "X", "Y"))
chr <- chr[! is.na(chr$CHROM), ]

pdf("./genome.pdf", width = 12, height = 20)  # density per chromosome

ggplot() + # geom_density() +
  geom_segment(data = chr, aes(x = 1, y = 1e2, xend = V2, yend = 1e2), size = 20, color = "lightblue") +
  geom_histogram(data = var, aes(x = POS), binwidth = 1e5, fill = "red") +
  theme_bw() + xlab("") + ylab("") + facet_grid(CHROM ~ .) # + coord_cartesian(ylim=c(1, 5e2))

dev.off()

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

CNA <- read.table("wes/BXSB_CNA.copynumber.called", header = T)
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

vcf_bxsb = vcf_bxsb[fixed(vcf_bxsb)$FILTER == "PASS"]
vcf_bxsb = vcf_bxsb[width(vcf_bxsb) == 1] # SNP only
vcf_bxsb = dropSeqlevels(vcf_bxsb, "chrM")
vcf_bxsb = renameSeqlevels(vcf_bxsb, gsub("chr", "", seqlevels(vcf_bxsb)))

idx <- sapply(split(mcols(all)$QUERYID, mcols(all)$GENEID), unique)

# summarize variant location by gene

variantByLocation <- sapply(names(idx), function(name) {
  d = all[mcols(all)$GENEID %in% name, c("QUERYID", "LOCATION")]
  table(mcols(d)$LOCATION[duplicated(d) == F])
})

rowSums(variantByLocation)

vep = ensemblVEP(file = sf)

sf = "/data/cgd/Sanger/REL-1505/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
(sanger_hdr = scanVcfHeader(sf))
ranges = GRanges(seqnames = 10, ranges = IRanges(start = 4455000, end = 4456000))
ranges = rowData(vcf_bxsb)
param = ScanVcfParam(geno = c("GT", "FI"), which = ranges)
vcf_sanger = readVcf(file = sf, genome = "mm10", param = param)

sanger_geno <- geno(vcf_sanger)
sanger_geno$GT
sanger_geno$FI
