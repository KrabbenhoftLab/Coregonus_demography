#!/bin/bash

#depth calculateions for new genome assembly mapping.

samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CA01.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CA01_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CA02.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CA02_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CA03.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CA03_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CA04.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CA04_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CH01.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CH01_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CH02.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CH02_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CH03.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CH03_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CH04.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CH04_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CK01.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CK01_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CK02.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CK02_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CK03.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CK03_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CK04.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CK04_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CN01.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CN01_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CN02.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CN02_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CZ01.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CZ01_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CZ02.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CZ02_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CZ03.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CZ03_depth.txt
samtools depth /mnt/krab1/NJCB_data/Coregonus_PSMC_v2/Cisco-illumina-bams/CZ04.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > CZ04_depth.txt