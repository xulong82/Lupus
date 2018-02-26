#!/bin/bash

dir="/home/xwang/Dropbox/GitHub/Lupus/WES"

file1=BXSB_GES13_04091_CGATGT_L005_trim
file2=nm3848_GES13_05080_ATCACG_L005_trim
file3=nm4690_GES13_04100_TAGCTT_L007_trim

qsub -v arg=$file1 $dir/bowtie.sh
qsub -v arg=$file2 $dir/bowtie.sh
qsub -v arg=$file3 $dir/bowtie.sh

