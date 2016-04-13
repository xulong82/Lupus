#!/bin/bash

module load Anaconda
source activate g2gtools

REF="/data/kbchoi/data/mm10/R75-REL1410/REF/C57BL6J.fa"
STRAIN="SB"
VCF_SNPS="/data/xwang/Lupus/emase/Sanger.UNC.Combined.SNPs.mm10.SB.vcf.gz"
GTF="/data/kbchoi/data/mm10/R75-REL1410/REF/C57BL6J.gtf"
DB="/data/kbchoi/data/mm10/R75-REL1410/REF/C57BL6J.gtf.db"

cd /data/xwang/Lupus/emase

# Patch SB genome from SNP

g2gtools patch -i ${REF} -s ${STRAIN} -v ${VCF_SNPS} -o ${STRAIN}.patched.fa

# Extract genes

g2gtools extract --genes -i ${STRAIN}.patched.fa -db ${DB} > ${STRAIN}.patched.genes.fa
g2gtools extract --genes -i ${REF} -db ${DB} > B6.genes.fa

g2gtools extract --transcripts -i ${STRAIN}.patched.fa -db ${DB} > ${STRAIN}.patched.transcripts.fa
g2gtools extract --transcripts -i ${REF} -db ${DB} > B6.transcripts.fa

