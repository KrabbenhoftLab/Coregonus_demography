##############################################################################################################################

# script to plot PSMC & SMC++ Coregonus
# author: NJCB
# njbacken@buffalo.edu

##############################################################################################################################

library(RColorBrewer)
library(rgdal)
library(gridExtra)
library(scales)
library(raster)
library(ggplot2)
library(grid)
library(dplyr)    
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(scales)
library(shades)


# specify the following parameters
wd <- "C:/Users/NJCB/Documents/Buffalo/Coregonus/Coregonus_PSMC-SMC_v2/results/" # top level working directory
setwd(wd)
getwd()
list.files()


########################################################
#                       PSMC                           #
########################################################
# set plotting values for PSMC

ls01_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS01_bs.txt"
ls01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS01_main.txt"
ls02_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS02_bs.txt"
ls02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS02_main.txt"
ls05_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS05_bs.txt"
ls05_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS05_main.txt"
ls07_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS07_bs.txt"
ls07_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS07_main.txt"
ls08_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS08_bs.txt"
ls08_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS08_main.txt"
ls09_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS09_bs.txt"
ls09_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS09_main.txt"


# read in PSMC data
ls01_p_bs <- read.delim(ls01_p_bs_file, header = TRUE, sep = "\t")
ls01_p_d <- read.delim(ls01_p_file, header = TRUE, sep = "\t")
ls02_p_bs <- read.delim(ls02_p_bs_file, header = TRUE, sep = "\t")
ls02_p_d <- read.delim(ls02_p_file, header = TRUE, sep = "\t")
ls05_p_bs <- read.delim(ls05_p_bs_file, header = TRUE, sep = "\t")
ls05_p_d <- read.delim(ls05_p_file, header = TRUE, sep = "\t")
ls07_p_bs <- read.delim(ls07_p_bs_file, header = TRUE, sep = "\t")
ls07_p_d <- read.delim(ls07_p_file, header = TRUE, sep = "\t")
ls08_p_bs <- read.delim(ls08_p_bs_file, header = TRUE, sep = "\t")
ls08_p_d <- read.delim(ls08_p_file, header = TRUE, sep = "\t")
ls09_p_bs <- read.delim(ls09_p_bs_file, header = TRUE, sep = "\t")
ls09_p_d <- read.delim(ls09_p_file, header = TRUE, sep = "\t")


# remove rows with "0" timepoints
ls01_p_bs <- ls01_p_bs[!(ls01_p_bs$x==0),] 
ls01_p_d <- ls01_p_d[!(ls01_p_d$x==0),] 
ls02_p_bs <- ls02_p_bs[!(ls02_p_bs$x==0),] 
ls02_p_d <- ls02_p_d[!(ls02_p_d$x==0),] 
ls05_p_bs <- ls05_p_bs[!(ls05_p_bs$x==0),] 
ls05_p_d <- ls05_p_d[!(ls05_p_d$x==0),] 
ls07_p_bs <- ls07_p_bs[!(ls07_p_bs$x==0),] 
ls07_p_d <- ls07_p_d[!(ls07_p_d$x==0),] 
ls08_p_bs <- ls08_p_bs[!(ls08_p_bs$x==0),] 
ls08_p_d <- ls08_p_d[!(ls08_p_d$x==0),] 
ls09_p_bs <- ls09_p_bs[!(ls09_p_bs$x==0),] 
ls09_p_d <- ls09_p_d[!(ls09_p_d$x==0),] 

#Scale y axis by multiplying by 10^4
ls01_p_bs$y <- ls01_p_bs$y*10^4 
ls01_p_d$y <- ls01_p_d$y*10^4 
ls02_p_bs$y <- ls02_p_bs$y*10^4 
ls02_p_d$y <- ls02_p_d$y*10^4 
ls05_p_bs$y <- ls05_p_bs$y*10^4 
ls05_p_d$y <- ls05_p_d$y*10^4 
ls07_p_bs$y <- ls07_p_bs$y*10^4 
ls07_p_d$y <- ls07_p_d$y*10^4 
ls08_p_bs$y <- ls08_p_bs$y*10^4 
ls08_p_d$y <- ls08_p_d$y*10^4 
ls09_p_bs$y <- ls09_p_bs$y*10^4 
ls09_p_d$y <- ls09_p_d$y*10^4 

########################################################
#                       SMC++                          #
########################################################

# set plotting values for SMC++
lean_s_file <- "./SMC++_MASKED/SN-lean_7.26e-9_t10-10000.gen.csv"
lean_s_bs_file <- "./SMC++_MASKED/SN-lean_bs_all.gen.csv"
sis_s_file <- "./SMC++_MASKED/SN-siscowet_7.26e-9_t10-10000.gen.csv"
sis_s_bs_file <- "./SMC++_MASKED/SN-siscowet_bs_all.gen.csv"


#read in data
lean_s_d <- read.csv(lean_s_file, header =TRUE)
lean_s_bs <- read.csv(lean_s_bs_file, header =TRUE)
sis_s_d <- read.csv(sis_s_file, header =TRUE)
sis_s_bs <- read.csv(sis_s_bs_file, header =TRUE)

#eliminate values with 0 time points
lean_s_d <- lean_s_d[!(lean_s_d$x==0),] 
lean_s_bs <- lean_s_bs[!(lean_s_bs$x==0),] 
sis_s_d <- sis_s_d[!(sis_s_d$x==0),] 
sis_s_bs <- sis_s_bs[!(sis_s_bs$x==0),] 

# scale x axis based on generation time
gen_time <- 16

lean_s_d$x_year <- lean_s_d$x * gen_time 
lean_s_bs$x_year <- lean_s_bs$x * gen_time 
sis_s_d$x_year <- sis_s_d$x * gen_time 
sis_s_bs$x_year <- sis_s_bs$x * gen_time 

# color blind safe palette
# color order matches alphabetical order of species
# i.e. CA, CH, CK, CN, CZ
colors <- c("#8b0000", "#8b0000", "#8b0000")
bsColor1 <- c("#EA5869") # make brighter version of the base color (good only for black)
bsColor2 <- c("#FAD000")
bsColor3 <- c("#000080")
bsColor4 <- c("#000000")
bsColor11 <- c("#8b0000") #*
bsColor14 <- c("#FF0000") #*
# geologic time boundaries, from ICS v2020/03 (https://stratigraphy.org/chart)
pliestocene <- c(11700, Inf) # use Inf to plot properly
# last glacial maximum, from Clark et al. 2009, https://science.sciencemag.org/content/325/5941/710.abstract
LGM <- c(19000, 33000)
# ice retreat from Great Lakes, from Dalton et al. 2020 https://www.sciencedirect.com/science/article/pii/S0277379119307619
deglaciation <- c(13500, 15500)



#plot PSMC

a <- ggplot(ls01_p_d, aes(x=x, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = ls01_p_d, aes(x=x, y=y, color=label), linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls01_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor11, alpha=0.1)) +
  geom_line(data = ls02_p_d, aes(x=x, y=y, color=label), linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls02_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor11, alpha=0.1)) +
  geom_line(data = ls07_p_d, aes(x=x, y=y, color=label), linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls07_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor11, alpha=0.1)) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  scale_color_manual(name="Species", values=colors) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL,
    limits=c(1000, 2800000)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL
  ) +
  annotation_logticks(sides="bl") + # log scale tick marks for both axes
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(3000,1000000))

a

c <- ggplot(ls01_p_d, aes(x=x, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = ls05_p_d, aes(x=x, y=y, color=label), color = bsColor14, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls05_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor14, alpha=0.1)) +
  geom_line(data = ls08_p_d, aes(x=x, y=y, color=label), color = bsColor14, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls08_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor14, alpha=0.1)) +
  geom_line(data = ls09_p_d, aes(x=x, y=y, color=label), color = bsColor14, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls09_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor14, alpha=0.1)) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  scale_color_manual(name="Species", values=colors) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL,
    limits=c(1000, 2800000)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL
  ) +
  annotation_logticks(sides="bl") + # log scale tick marks for both axes
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(3000,1000000))

c



#plot SMC++

b <- ggplot(lean_s_d, aes(x=x_year, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(linewidth = 1.25, color = bsColor11) +
  geom_line(data = lean_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor11, alpha=0.2)) +
  geom_line(data = lean_s_d, aes(x=x_year, y=y, group=label), color=bsColor11) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  scale_color_manual(name="Species", values=colors) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL,
    limits=c(1000, 2800000)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL
  ) +
  annotation_logticks(sides="bl") + # log scale tick marks for both axes
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(3000,1000000))

b

d <- ggplot(sis_s_d, aes(x=x_year, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(linewidth = 1.25, color = bsColor14) +
  geom_line(data = sis_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor14, alpha=0.2)) +
  geom_line(data = sis_s_d, aes(x=x_year, y=y, group=label), color=bsColor14) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  scale_color_manual(name="Species", values=colors) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL,
    limits=c(1000, 2800000)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL
  ) +
  annotation_logticks(sides="bl") + # log scale tick marks for both axes
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(3000,1000000))

d
  


plots <- plot_grid(a,
                   b,
                   c,
                   d,
                   ncol = 2)

  
plots  
  

wdII <- ("C:/Users/NJCB/Documents/Buffalo/Coregonus/Coregonus_PSMC-SMC_v2/figures_and_tables/SupFig5_Snamaycush_PSMC_SMC_bs_plot_grid/")
ggsave(paste(wdII, "./", "SupFig5_Snamaycush_PSMC_SMC_bs_plot_grid_11-21-23", ".png", sep=""), plot=plots, units = "in", dpi=1000, height=8, width=10)
ggsave(paste(wdII, "./", "SupFig5_Snamaycush_PSMC_SMC_bs_plot_grid_11-21-23", ".svg", sep=""), plot=plots, units = "in", dpi=1000, height=8, width=10)
