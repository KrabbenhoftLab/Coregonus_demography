#!/bin/bash

HOME="/mnt/krab1/NJCB_data/Coregonus_PSMC_v2/bam_stats"
cd ${HOME}

FILELIST=`cat temp`
for f in $FILELIST
do
	samtools flagstat -@ 8 /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/${f}.bam > ${HOME}/${f}_flagstat.txt
	samtools mpileup /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/${f}.bam | awk -v X="${1}" '$4>=X' | wc -l > ${HOME}/${f}_breadth.txt
done
