#!/bin/bash
#PBS -l mem=128gb,nodes=1:ppn=20,walltime=10:00:00

module load bwa

cd /data/xwang
bwa index mm10.fa

bwa aln -t 20 ref.fa -b1 reads.bam > 1.sai 
bwa aln -t 20 ref.fa -b2 reads.bam > 2.sai 
bwa sampe ref.fa 1.sai 2.sai reads.bam reads.bam > aln.sam

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
       ../bam/"$file".sam
 
