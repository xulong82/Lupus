#!/bin/bash

dir1="/data/xwang/Lupus/WES/TRIM"
dir2="/home/xwang/Dropbox/GitHub/Lupus/WES"

file1=BXSB_GES13_04091_CGATGT
file2=nm3848_GES13_05080_ATCACG

qsub -v arg=$file1 $dir2/bowtie.sh
qsub -v arg=$file2 $dir2/bowtie.sh

# files=`find $dir1 -name '*X_R1.fastq'`
# 
# for name1 in $files; do
#   name2=`basename $name1`
#   name3=${name2/_R1.fastq/}
#   echo $name3
#   qsub -v arg=$name3 $dir2/bowtie.sh
# done
  
