#!/bin/sh
#AISO Script (ABOUT, INPUT, STEPS, OUTPUT) for running psmc analysis.
#TJ Krabbenhoft, May 2020, University at Buffalo
#This version run on MacPro

#ABOUT:
#This script takes input bam files from a genome assembly and runs PSMC analysis. The analysis requires a diploid consensus (i.e., SNPs called throughout the genome, with files created in Step 2).

#INPUT:
SEQDIR="/mnt/krab1/NJCB_data/Snamaycush_PSMC_v2/original_VCFs"

cd ${SEQDIR}
FILELIST=`cat temp`

#STEP 2: Consenus calling:
for f in ${FILELIST}
do
#set variables for math
depth=`cat /mnt/krab1/NJCB_data/Snamaycush_PSMC_v2/avg_depth/${f}_depth.txt`
min=`printf %.0f $(echo "$depth / 3" | bc -l)`
max=`printf %.0f $(echo "$depth * 2" | bc -l)`
vcfutils.pl vcf2fq -d ${min} -D ${max} ${f}_temp1.vcf | gzip > ${f}.diploid.fq.gz
echo "Concensus calling completed. It won't be long now."
done
