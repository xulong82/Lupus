#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=1,walltime=10:00:00

module load bowtie
module load samtools
module load Anaconda

source activate emase

# file=${arg}

# cd /data/xwang/Lupus/emase
# create-hybrid -F B6.genes.fa,SB.patched.genes.fa -s B,S --create-bowtie-index
 
# hybrid genome
REF_H="/data/xwang/Lupus/emase/emase.pooled.targets.bowtie1"

# Info for the B6 reference (ultimate reference)
INFO1="/data/xwang/Lupus/emase/emase.transcripts.info"
# Info for the hybrid references 
INFO2="/data/xwang/Lupus/emase/emase.pooled.targets.info"

FASTQ="/data/xwang/Lupus/wes/fastq/BXSB_GES13_04091_CGATGT_L005_trim_R1.fastq"
DIREC="/data/xwang/Lupus/wes/fastq"
BASE1=`basename $FASTQ`
BASE2=${BASE1/_R1.fastq/}
 
# alignment in parallel (ppn=20)

cd /data/xwang/Lupus/wes/emase
bowtie -p 20 -q -a --best --strata --sam -v 3 ${REF_H} ${DIREC}/${BASE2}_R1.fastq ${BASE2}_R1.sam
bowtie -p 20 -q -a --best --strata --sam -v 3 ${REF_H} ${DIREC}/${BASE2}_R2.fastq ${BASE2}_R2.sam

# sam to bam format

samtools view -bS ${BASE2}_R1.sam > ${BASE2}_R1.bam
samtools view -bS ${BASE2}_R2.sam > ${BASE2}_R2.bam
  
# bam to emase format

bam-to-emase -a ${BASE2}_R1.bam -i ${INFO1} -s B,S -o ${BASE2}_R1.h5
bam-to-emase -a ${BASE2}_R2.bam -i ${INFO1} -s B,S -o ${BASE2}_R2.h5
  
# join paired-files

combine-emase-files -i ${BASE2}_R1.h5,${BASE2}_R2.h5 -o ${BASE2}.h5
  
# emase estimation of expression level

run-emase -i ${BASE2}.h5 -L ${INFO2} -M 4 -o emase -c

