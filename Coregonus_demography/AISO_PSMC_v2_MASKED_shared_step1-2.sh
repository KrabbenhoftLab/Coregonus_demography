#!/bin/sh
#AISO Script (ABOUT, INPUT, STEPS, OUTPUT) for running psmc analysis.
#TJ Krabbenhoft, May 2020, University at Buffalo


#ABOUT:
#This script takes input bam files from a genome assembly and runs PSMC analysis. The analysis requires a diploid consensus (i.e., SNPs called throughout the genome, with files created in Step 2).

#INPUT:

HOME="/mnt/sdc1/Coregonus_PSMC-SMC/PSMC_masked_genome_04-10-23/"
SEQDIR="/mnt/sdc1/Coregonus_PSMC-SMC/Cisco-illumina-genomes/"
BWADIR=""
REF="/mnt/sdc1/Coregonus_PSMC-SMC/PSMC_masked_genome_04-10-23/CA_GEN_1-flye_III-pilon2-purged.FINAL.chromonomer-03-21-2022.masked.fasta" #Reference Assembly
PSMCDIR="/home/krablab/Documents/apps/psmc-master" #PSMC directory


#STEP 1: READ MAPPING
cd ${HOME}
${BWADIR}bwa index ${REF}
FILELIST=`cat temp`
#FILELIST="CZ03" #This line just for running individual samples; otherwise use line above.
for f in ${FILELIST}
do
${BWADIR}bwa mem -t 8 ${REF} ${SEQDIR}${f}_1.fq.gz ${SEQDIR}${f}_2.fq.gz | samtools sort -@ 8 -o ${f}.bam -O bam
samtools index -@ 8 ${f}.bam
done


#STEP 2: Consenus calling:

for f in $FILELIST
do
bcftools mpileup -C50 --threads 40 -f ${REF} ${f}.bam > temp.vcf
bcftools call --threads 40 -c temp.vcf > ${f}_temp1.vcf
rm temp.vcf
bedtools intersect -v -a ./${f}_temp1.vcf.gz -b ../CA_GEN_1_assembly/CA_GEN_1.final_repeat_mask.complex.reformat_4_column.filtered.tabbed.bed -header > /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/masked_vcf/${f}_temp1.masked.vcf
done


cd /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/masked_vcf

for f in ${FILELIST}
do
#echo "SNP calling completed! Things are looking up!"
#set variables for math
depth=`cat /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/avg_depth/${f}_depth.txt`
min=`printf %.0f $(echo "$depth / 3" | bc -l)`
max=`printf %.0f $(echo "$depth * 2" | bc -l)`
vcfutils.pl vcf2fq -d ${min} -D ${max} ${f}_temp1.masked.vcf | gzip > ${f}.masked.diploid.fq.gz
echo "Consensus calling completed."
done
