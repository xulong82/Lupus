#!/bin/bash
#PBS -l mem=128gb,nodes=1:ppn=20,walltime=10:00:00

module load bowtie
module load samtools

# C57BL6J WES: Laura Reinholdt
# /data/gananda/projects/ExomeSeq_Reinholdt/inbred_strains/
# Not the same library; not the same sequencer

# Heather Fairfield's sample matches the most
# /data/hef/data_from_laura

# echo "Build bowtie index for UCSC mm10 whole genome."
# cd /data/xwang/REF
# bowtie-build Mm10Genome.fa Mm10Genome

cd /data/xwang/Lupus/WES/TRIM
ref="/data/xwang/REF/Mm10Genome"
file=${arg}
  
bowtie -p 20 \
       -q \
       --phred33-quals \
       --chunkmbs 512 \
       --sam \
       -a --best --strata \
       "$ref" \
       -1 "$file"_R1.fastq \
       -2 "$file"_R2.fastq \
       "$file".sam
 
