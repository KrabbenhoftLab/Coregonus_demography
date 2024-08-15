## This repository contains scripts for analyses presented in Backenstose et al., 2024. 
Title: *"Origin of the Laurentian Great Lakes fish fauna through upward adaptive radiation cascade prior to the Last Glacial Maximum"*

Related data files can be found here: https://doi.org/10.5061/dryad.n02v6wx59 and here https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1062807

## Directory description
### Coregonus_IQTree2_Fig1
-Script for used for alignment of variable sites within Coregonus and generation of ML tree.

### Coregonus_PCA_Fig1_SuppFig1
-Script and sample data to perform PCA analyses using genome-wide SNP data.

### Coregonus_demography
**Alignment to reference and SNP calling**
  AISO_PSMC_v2_MASKED_shared_step1-2.sh
**PSMC analyses with bootstrapping**
  AISO_PSMC_v2_MASKED_shared-step3_param7_bootstrap.sh
**Post processing of PSMC output, with mutation rate of 7.26e-09 mutations/site/generation and a generation time of 6 years**
  all_species_PSMC_MASKED_shared_bootstrap_plotting_11-13-2023.sh
**Post processing of PSMC output, with variable mutation rates**
  all_species_PSMC_MASKED_MUT_RATE_shared_bootstrap_plotting_11-30-2023.sh
**Alignment statistics
  bam_stats.sh
  calc_average_depth.sh


## Please cite:
Backenstose, N.J.C., MacGuigan, D.J., Osborne, C.A. et al. Origin of the Laurentian Great Lakes fish fauna through upward adaptive radiation cascade prior to the Last Glacial Maximum. Commun Biol 7, 978 (2024). https://doi.org/10.1038/s42003-024-06503-z

