#!/bin/sh
#PBS -l mem=128gb,nodes=1:ppn=20,walltime=10:00:00

picard="$HOME/picard-tools-1.139/picard.jar"

cd /data/xwang/Lupus/WES/bam

fastq="BXSB_GES13_04091_CGATGT_L005_trim"
bamfile="BXSB_GES13_04091_CGATGT_L005_trim.sorted.bam"

# cd /data/xwang/Lupus/WES/fastq
# 
# java -jar $picard FastqToSam \
#           FASTQ=${fastq}_R1.fastq FASTQ2=${fastq}_R2.fastq OUTPUT=${fastq}.ubam \
# 	  SAMPLE_NAME=BXSB
# 
# java -jar $picard AddOrReplaceReadGroups \
#           INPUT=$bamfile OUTPUT=${bamfile/.bam/_rg.bam/} \
# 	  RGID=1 RGLB=1 RGPL=ILLUMINA RGSM=BXSB RGPU=1
# 
  java -jar $picard FixMateInformation \
            INPUT=${bamfile/.bam/_rg.bam/} OUTPUT=${bamfile/.bam/_rg_fix.bam/}
# 
# java -jar $picard ValidateSamFile INPUT=${bamfile/.bam/_rg_fix.bam/}
# 
# # mark duplicates
# 
# java -jar $picard MarkDuplicates.jar \ 
#     	  INPUT=${bamfile/.bam/_rg.bam/} \
#           OUTPUT=${bamfile/.bam/_rg_dedup.bam/} \
# 	  METRICS_FILE=${bamfile/.bam/_rg_dedup.txt/}
# 
# java -jar BuildBamIndex.jar \ 
#     INPUT=dedup_reads.bam 
# 
# echo $0
# 
# module load java/1.7.0
# 
# dir="/data/xwang/Lupus/WES/fastq"
# 
# files=`find $dir -name '*_R1_ALL.fastq'`
# 
# for name1 in $files; do
#   name2=`basename $name1`
#   name3=${name2/_R1_ALL.fastq/}
#   echo $name3
#   java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 20 -phred33 \
#             "$dir"/"$name3"_R1_ALL.fastq \
#             "$dir"/"$name3"_R2_ALL.fastq \
#             "$dir"/"$name3"_trim_R1.fastq \
#             "$dir"/"$name3"_unpaired_R1.fastq \
#             "$dir"/"$name3"_trim_R2.fastq \
#             "$dir"/"$name3"_unpaired_R2.fastq \
#             LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:60
# done
# 
