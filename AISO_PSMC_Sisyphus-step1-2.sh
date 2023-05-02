#!/bin/sh
#AISO Script (ABOUT, INPUT, STEPS, OUTPUT) for running psmc analysis. 
#TJ Krabbenhoft, May 2020, University at Buffalo
#This version run on MacPro

#ABOUT:
#This script takes input bam files from a genome assembly and runs PSMC analysis. The analysis requires a diploid consensus (i.e., SNPs called throughout the genome, with files created in Step 2). 

#INPUT:
SEQDIR="/mnt/sdc1/Coregonus_PSMC-SMC/Cisco-illumina-genomes/"
BAMSUF=".bwa.sorted_pdrm.bam"  #Mapped reads bam file suffix. 
BWADIR=""
REF="/mnt/sdc1/Coregonus_PSMC-SMC/CA_GEN_1-shasta_1kb-pilon3_chromonomer_filtered.fasta" #Reference Assembly
PSMCDIR="/home/krablab/Documents/apps/psmc-master" #PSMC directory


#STEP 1: READ MAPPING
cd ${SEQDIR}
${BWADIR}bwa index ${REF}
FILELIST=`cat temp1`
#FILELIST="CZ03" #This line just for running individual samples; otherwise use line above.
for f in ${FILELIST}
do
${BWADIR}bwa mem -t8 ${REF} ${f}_1.fq.gz ${f}_2.fq.gz | samtools sort -o ${f}.bam -O bam
samtools index ${f}.bam
done

#STEP 2: Consenus calling:
for f in ${FILELIST}
do
bcftools mpileup -C50 -f ${REF} ${f}.bam > temp.vcf
echo "Read mapping has finished. Boy am I glad that is over."
bcftools call -c temp.vcf > temp1.vcf
#rm temp.vcf
echo "SNP calling completed! Things are looking up!"
vcfutils.pl vcf2fq -d 10 -D 120 temp1.vcf | gzip > ${f}.diploid.fq.gz
#rm temp1.vcf
echo "Concensus calling completed. It won't be long now."
done
