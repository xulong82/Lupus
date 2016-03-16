#!/bin/bash

# Make VCF manually

input="/data/vmp/snps/Sanger.UNC.Combined.SNPs.mm10.SB.txt"
input="/data/vmp/snps/Sanger.UNC.Combined.SNPs.mm10.BXSB.txt"

file=`basename $input`

cd /data/xwang/Lupus/emase

echo "Adjust in the vcf order"
awk 'BEGIN { FS="\t" } { print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6 }' $input > temp1

awk -v OFS="\t" 'BEGIN { FS="\t" } { if($3!~/rs/) $3="."; print $0 }' temp1 > temp2

echo "Remove header"
sed '1d' temp2 > temp3

echo "Extra columns"
awk -v OFS="\t" 'BEGIN { FS="\t" } { $7="PASS"; $8="N"; $9="GT:FI"; $10="1/1:1" } { print $0 }' temp3 > temp4

echo "VCF header"
awk 'BEGIN { print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSB" } { print $0 }' temp4 > temp5

echo "Keep only mutated"
awk 'BEGIN { FS="\t" } { if ($4 != $5) print $0 }' temp5 > temp6

echo "Remove the non-called"
awk 'BEGIN { FS="\t" } { if ($5 != "N") print $0 }' temp6 > temp7

cp temp7 ${file/.txt/.vcf}

rm temp[1-7]

bgzip ${file/.txt/.vcf}
tabix ${file/.txt/.vcf.gz}

