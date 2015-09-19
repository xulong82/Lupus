#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=20,walltime=10:00:00

file=${arg}

module load rsem
module load bowtie/1.0.0

# BUILD C3H TRANSCRIPTOME
# cd /data/xwang/Testis/REF
# rsem-prepare-reference C3H.transcripts_ERCC.fa C3H.transcripts
# rsem-prepare-reference beth.fa beth.transcripts 

dir1="/data/xwang/Lupus/BXSB/trim"

# ref="/data/xwang/Lupus/Anuj/RSEM_ref/BXSB_pseudo"
# dir2="/data/xwang/Lupus/RSEM_BXSB"

dir2="/data/xwang/Lupus/RSEM_B6"
ref="/data/xwang/RSEM/GRCm38"

rsem-calculate-expression -p 20 \
			  --bowtie-phred33-quals \
   			  --forward-prob 0.5 \
   			  --paired-end \
                          "$dir1"/"$file"_R1.fastq \
                          "$dir1"/"$file"_R2.fastq \
   			  "$ref" \
   			  "$dir2"/"$file"

