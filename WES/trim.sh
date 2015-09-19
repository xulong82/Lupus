#!/bin/sh
#PBS -l mem=128gb,nodes=1:ppn=20,walltime=10:00:00

# Read trimming with Trimmomatic

echo $0

module load java/1.7.0

dir1="/data/xwang/Lupus/WES/RAW"
dir2="/data/xwang/Lupus/WES/TRIM"

files=`find $dir1 -name '*_R1_ALL.fastq'`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/_L005_R1_ALL.fastq/}
  echo $name3
  java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 20 -phred33 \
            "$dir1"/"$name3"_L005_R1_ALL.fastq \
            "$dir1"/"$name3"_L005_R2_ALL.fastq \
            "$dir2"/"$name3"_R1.fastq \
            "$dir2"/"$name3"_unpaired_R1.fastq \
            "$dir2"/"$name3"_R2.fastq \
            "$dir2"/"$name3"_unpaired_R2.fastq \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:60
done

