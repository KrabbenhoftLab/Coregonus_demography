#!/bin/bash

# get alignment containing only variable sites
iqtree2 --seqtype DNA -s Coregonus_PCA_biallelic_no_missing.maf0.05.noSingletons.min4.phy -m MFP+ASC -B 1000 -bnni -T AUTO --threads-max 8 --prefix Coregonus_PCA_biallelic_no_missing.maf0.05.noSingletons.min4 > iqtree.out 2> iqtree.err

# run IQTree2 using only variable sites
iqtree2 --seqtype DNA -s Coregonus_PCA_biallelic_no_missing.maf0.05.noSingletons.min4.varsites.phy -m MFP+ASC -B 1000 -bnni -T 12 --prefix Coregonus_PCA_biallelic_no_missing.maf0.05.noSingletons.min4.varsites > iqtree.2.out 2> iqtree.2.err
