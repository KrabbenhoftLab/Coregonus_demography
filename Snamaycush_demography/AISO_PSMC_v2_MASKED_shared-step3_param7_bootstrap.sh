#!/bin/sh
#AISO Script (ABOUT, INPUT, STEPS, OUTPUT) for running psmc analysis. 
#TJ Krabbenhoft, May 2020, University at Buffalo

#ABOUT:
#This script takes input bam files from a genome assembly and runs PSMC analysis. The analysis requires a diploid consensus (i.e., SNPs called throughout the genome, with files created in Step 2). 

#INPUT:
SEQDIR="/mnt/krab1/NJCB_data/Snamaycush_PSMC_v2/masked_vcf"
PSMCDIR="/home/krablab/Documents/apps/psmc" #PSMC directory

cd ${SEQDIR}
FILELIST=`cat temp`

OUT="/mnt/krab1/NJCB_data/Snamaycush_PSMC_v2/PSMC_MASKED_param7_bootstrap_raw_bam"
mkdir -p ${OUT}
cd ${OUT}
#STEP 3: Run PSMC:
for f in ${FILELIST}
do
${PSMCDIR}/utils/fq2psmcfa -q20 ${SEQDIR}/${f}.masked.diploid.fq.gz > ${f}.single.diploid.psmcfa
${PSMCDIR}/psmc -N25 -t5 -r5 -p "4+20*2+6*4+4" -o ${f}.single.diploid.psmc ${f}.single.diploid.psmcfa
${PSMCDIR}/utils/splitfa ${f}.single.diploid.psmcfa > ${f}.split.psmcfa
seq 100 | xargs -t -n 1 -P 50 -i ${PSMCDIR}/psmc -N25 -t5 -r5 -b -p "4+20*2+6*4+4" -o ${f}_round-{}.psmc ${f}.split.psmcfa | sh
cat ${f}.single.diploid.psmc ${f}_round-*.psmc > ${f}.combined.psmc
done

echo "bootstrapping complete"

#OUTPUT: 
#Relevant output files include: 
#Image file with PSMC graph
#Make PSMC plot with multiple samples and correct for missing data.

