#!/bin/bash
#PBS -l mem=128gb,nodes=1:ppn=20,walltime=10:00:00

module load perl
path_vep="/home/xwang/ensembl-tools-release-81/scripts/variant_effect_predictor"
cd $HOME/Dropbox/GitHub/Lupus/WES 

input=BXSB_SNP_filtered
perl $path_vep/variant_effect_predictor.pl --fork 20 \
     --gmaf --symbol --protein --biotype --regulatory --force_overwrite \
     -i ./$input.vcf -o ./$input_vep.txt

