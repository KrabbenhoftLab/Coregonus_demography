#!/bin/bash

HOME="/mnt/krab1/NJCB_data/Coregonus_PSMC_v2/vcfs_PCA/"

cd ${HOME}

bcftools merge --threads 15 -Oz -l vcf_list.txt > Coregonus_PCA.vcf.gz

bcftools view -e 'GT[*]="mis"' --threads 15 -m 2 -M 2 -Oz Coregonus_PCA.vcf.gz > Coregonus_PCA_biallelic_no_missing.vcf.gz


