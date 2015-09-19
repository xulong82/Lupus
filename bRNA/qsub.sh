#!/bin/bash

dir1="/data/xwang/Lupus/BXSB/trim"
dir2="/home/xwang/Dropbox/GitHub/Lupus/bRNA"

files=`find $dir1 -name '*L005_R1.fastq'`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/_R1.fastq/}
  echo $name3
  qsub -v arg=$name3 $dir2/rsem.sh
done
  
# files=`find $dir2 -name '*.sh.e*'`
# for name1 in $files; do
#   name2=`basename $name1`
#   echo $name2 >> 1
#   grep "at least" $name2 >> 1
# done 

