#!/bin/bash

dir1="/data/xwang/Lupus/bxsb_rna/trim"
dir2="/home/xwang/Dropbox/GitHub/Lupus/shell"

files=`find $dir1 -name '*_L005_R1.fastq'`

for name1 in $files; do
  name2=`basename $name1`
  echo $name2
  qsub -v arg=$name2 $dir2/emase_rna.sh
done
  
