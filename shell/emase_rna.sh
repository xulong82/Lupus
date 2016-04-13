#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=10,walltime=10:00:00

######################

module load bowtie
module load samtools
module load Anaconda

source activate emase

# cd /data/xwang/Lupus/emase
# create-hybrid -F B6.genes.fa,SB.patched.genes.fa -s B,S --create-bowtie-index
# create-hybrid -F B6.transcripts.fa,SB.patched.transcripts.fa -s B,S --create-bowtie-index
 
# hybrid genome
REF_H="/data/xwang/Lupus/emase/emase.pooled.targets.bowtie1"

# Info for the B6 reference (ultimate reference)
INFO1="/data/xwang/Lupus/emase/emase.transcripts.info"
# Info for the hybrid references 
INFO2="/data/xwang/Lupus/emase/emase.pooled.targets.info"

######################

FASTQ=${arg}
DIREC="/data/xwang/Lupus/bxsb_rna/trim"
BASEN=${FASTQ/_R1.fastq/}
 
cd /data/xwang/Lupus/bxsb_rna/emase

# alignment in parallel (ppn=10)
bowtie -p 10 -q -a --best --strata --sam -v 3 ${REF_H} ${DIREC}/${BASEN}_R1.fastq ${BASEN}_R1.sam
bowtie -p 10 -q -a --best --strata --sam -v 3 ${REF_H} ${DIREC}/${BASEN}_R2.fastq ${BASEN}_R2.sam

# sam to bam format
samtools view -bS ${BASEN}_R1.sam > ${BASEN}_R1.bam
samtools view -bS ${BASEN}_R2.sam > ${BASEN}_R2.bam
  
# bam to emase format
bam-to-emase -a ${BASEN}_R1.bam -i ${INFO1} -s B,S -o ${BASEN}_R1.h5
bam-to-emase -a ${BASEN}_R2.bam -i ${INFO1} -s B,S -o ${BASEN}_R2.h5
  
# join paired-files
combine-emase-files -i ${BASEN}_R1.h5,${BASEN}_R2.h5 -o ${BASEN}.h5
  
# emase estimation of expression level
run-emase -i ${BASEN}.h5 -L ${INFO2} -M 4 -o ${BASEN}.emase -c

