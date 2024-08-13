#!/bin/sh
#AISO Script (ABOUT, INPUT, STEPS, OUTPUT) for running psmc analysis. 
#TJ Krabbenhoft, May 2020, University at Buffalo
#This version run on MacPro

#ABOUT:
#This script takes input bam files from a genome assembly and runs PSMC analysis. The analysis requires a diploid consensus (i.e., SNPs called throughout the genome, with files created in Step 2). 

#INPUT:
SEQDIR="/mnt/krab1/NJCB_data/Coregonus_PSMC_v2/masked_vcf"
PSMCDIR="/home/krablab/Documents/apps/psmc" #PSMC directory

cd ${SEQDIR}
FILELIST=`cat temp`

OUT="/mnt/krab1/NJCB_data/Coregonus_PSMC_v2/PSMC_masked_param7_bootstrap_raw_bam"
mkdir -p ${OUT}
cd ${OUT}
#STEP 3: Run PSMC:
for f in ${FILELIST}
do
${PSMCDIR}/utils/fq2psmcfa -q20 ${SEQDIR}/${f}.masked.diploid.fq.gz > ${f}.single.diploid.psmcfa
${PSMCDIR}/psmc -N25 -t5 -r5 -p "4+20*2+6*4+4" -o ${f}.single.diploid.psmc ${f}.single.diploid.psmcfa
${PSMCDIR}/utils/splitfa ${f}.single.diploid.psmcfa > ${f}.split.psmcfa
seq 100 | xargs -t -n 1 -P 500 -i ${PSMCDIR}/psmc -N25 -t5 -r5 -b -p "4+20*2+6*4+4" -o ${f}_round-{}.psmc ${f}.split.psmcfa | sh
cat ${f}.single.diploid.psmc ${f}_round-*.psmc > ${f}.combined.psmc
done

echo "bootstrapping complete"
echo "It's finished. Go make some plots"

#OUTPUT: 
#Relevant output files include: 
#Image file with PSMC graph
#Make PSMC plot with multiple samples and correct for missing data:
#/Users/krablab/Downloads/psmc-master/utils/psmc_plot.pl -M "${f}=0.3,CA03=0.3,CH01=0.3,CH04=0.3,CK01=0.3,CK02=0.3,CZ01=0.3,CZ02=0.3,CN01=0.3,CN02=0.3" prefix ${f}.psmc CA03.psmc CH01.psmc CH04.psmc CK01.psmc CK02.psmc CZ01.psmc CZ02.psmc CN01.psmc CN02.psmc

