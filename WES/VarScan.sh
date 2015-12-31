#!/bin/bash
#PBS -l mem=128gb,nodes=1:ppn=1,walltime=10:00:00

module load bowtie
module load samtools
module load java

cd /data/xwang/Lupus/WES/bam

sams=`find ./ -name '*.sam'`
bams=`find ./ -name '*.sorted.bam'`
mm10="/data/xwang/REF/Mm10Genome"

# for idx in $sams; do
#   idx2=`basename $idx`
#   file=${idx2/.sam/}
#   echo $file
#   echo "sam to bam, sorting, indexing"
#   samtools view -bSF 4 -o ${file}.bam ${file}.sam
#   samtools sort ${file}.bam ${file}.sorted
#   samtools index ${file}.sorted.bam 
# done

cd /data/xwang/Lupus/WES/var

for idx in $bams; do
  idx2=`basename $idx`
  file=${idx2/.sorted.bam/}
  echo $file
  echo "coverage calling"
  samtools mpileup -q 1 -f ${mm10}.fa ../bam/${file}.sorted.bam > ${file}.mpileup
  awk '{if($4 >= 6) print $0}' ${file}.mpileup > ${file}_filtered.mpileup
done

# case="BXSB_GES13_04091_CGATGT_filtered.mpileup"
# control="nm3848_GES13_05080_ATCACG_filtered.mpileup"

# SNP
# java -Xmx128G -jar $HOME/VarScan.v2.3.9.jar pileup2snp $case --p-value 0.01 > BXSB_SNP.vsf
# java -Xmx128G -jar $HOME/VarScan.v2.3.9.jar mpileup2snp $case --output-vcf 1 --p-value 0.01 > BXSB_SNP.vcf
# java -Xmx128G -jar $HOME/VarScan.v2.3.9.jar pileup2snp $control --p-value 0.01 > nm3848_SNP.vsf
# INDEL
# java -Xmx128G -jar $HOME/VarScan.v2.3.9.jar pileup2indel $case --p-value 0.01 > BXSB_IND.vsf
# java -Xmx128G -jar $HOME/VarScan.v2.3.9.jar pileup2indel $control --p-value 0.01 > nm3848_IND.vsf

# java -jar $HOME/VarScan.v2.3.9.jar filter BXSB_SNP.vcf --min-reads2 50 > BXSB_SNP_filtered.vcf

# Copy number alteration
# awk '{(tot+=$4)}; END{print "total base:" tot}' $case
# awk '{(tot+=$4)}; END{print "total base:" tot}' $control

# java -Xmx128G -jar $HOME/VarScan.v2.3.9.jar copynumber $control $case BXSB_CNA --data-ratio 1.473508
# java -Xmx128G -jar $HOME/VarScan.v2.3.9.jar copyCaller BXSB_CNA.copynumber --output-file BXSB_CNA.copynumber.called


# NCBI Genome Reference Consortium (pseudo-autosome region)
# X: 169,969,759 -> 170,931,299
# Y: 90,745,845 -> 91,644,698
# samtools mpileup -r chrX:165000000-175000000 -f "$ref".fa "$file".sorted.bam > "$file"_chrX.pileup
# samtools mpileup -r chrY:85000000-95000000 -f "$ref".fa "$file".sorted.bam > "$file"_chrY.pileup
# awk '{print $1"\t"$2"\t"$3"\t"$4}' "$file"_chrX.pileup > chrX.pileup
# awk '{print $1"\t"$2"\t"$3"\t"$4}' "$file"_chrY.pileup > chrY.pileup

