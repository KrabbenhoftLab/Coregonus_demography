##############################################################################################################################

# script to plot PSMC & SMC++ Coregonus and Snamaycush
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
ca01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA01_main.txt"
ca02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA02_main.txt"
ca03_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA03_main.txt"
ca04_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA04_main.txt"
ch01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH01_main.txt"
ch02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH02_main.txt"
ch03_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH03_main.txt"
ch04_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH04_main.txt"
ck01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK01_main.txt"
ck02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK02_main.txt"
ck03_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK03_main.txt"
ck04_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK04_main.txt"
cn01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CN01_main.txt"
cn02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CN02_main.txt"
ls01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS01_main.txt"
ls02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS02_main.txt"
ls05_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS05_main.txt"
ls07_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS07_main.txt"
ls08_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS08_main.txt"
ls09_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/LS09_main.txt"


# read in PSMC data
ca01_p_d <- read.delim(ca01_p_file, header = TRUE, sep = "\t")
ca02_p_d <- read.delim(ca02_p_file, header = TRUE, sep = "\t")
ca03_p_d <- read.delim(ca03_p_file, header = TRUE, sep = "\t")
ca04_p_d <- read.delim(ca04_p_file, header = TRUE, sep = "\t")
ch01_p_d <- read.delim(ch01_p_file, header = TRUE, sep = "\t")
ch02_p_d <- read.delim(ch02_p_file, header = TRUE, sep = "\t")
ch03_p_d <- read.delim(ch03_p_file, header = TRUE, sep = "\t")
ch04_p_d <- read.delim(ch04_p_file, header = TRUE, sep = "\t")
ck01_p_d <- read.delim(ck01_p_file, header = TRUE, sep = "\t")
ck02_p_d <- read.delim(ck02_p_file, header = TRUE, sep = "\t")
ck03_p_d <- read.delim(ck03_p_file, header = TRUE, sep = "\t")
ck04_p_d <- read.delim(ck04_p_file, header = TRUE, sep = "\t")
cn01_p_d <- read.delim(cn01_p_file, header = TRUE, sep = "\t")
cn02_p_d <- read.delim(cn02_p_file, header = TRUE, sep = "\t")
ls01_p_d <- read.delim(ls01_p_file, header = TRUE, sep = "\t")
ls02_p_d <- read.delim(ls02_p_file, header = TRUE, sep = "\t")
ls05_p_d <- read.delim(ls05_p_file, header = TRUE, sep = "\t")
ls07_p_d <- read.delim(ls07_p_file, header = TRUE, sep = "\t")
ls08_p_d <- read.delim(ls08_p_file, header = TRUE, sep = "\t")
ls09_p_d <- read.delim(ls09_p_file, header = TRUE, sep = "\t")

# remove rows with "0" timepoints
ca01_p_d <- ca01_p_d[!(ca01_p_d$x==0),] 
ca02_p_d <- ca02_p_d[!(ca02_p_d$x==0),] 
ca03_p_d <- ca03_p_d[!(ca03_p_d$x==0),] 
ca04_p_d <- ca04_p_d[!(ca04_p_d$x==0),] 
ch01_p_d <- ch01_p_d[!(ch01_p_d$x==0),] 
ch02_p_d <- ch02_p_d[!(ch02_p_d$x==0),] 
ch03_p_d <- ch03_p_d[!(ch03_p_d$x==0),] 
ch04_p_d <- ch04_p_d[!(ch04_p_d$x==0),] 
ck01_p_d <- ck01_p_d[!(ck01_p_d$x==0),] 
ck02_p_d <- ck02_p_d[!(ck02_p_d$x==0),] 
ck03_p_d <- ck03_p_d[!(ck03_p_d$x==0),] 
ck04_p_d <- ck04_p_d[!(ck04_p_d$x==0),] 
cn01_p_d <- cn01_p_d[!(cn01_p_d$x==0),] 
cn02_p_d <- cn02_p_d[!(cn02_p_d$x==0),] 
ls01_p_d <- ls01_p_d[!(ls01_p_d$x==0),] 
ls02_p_d <- ls02_p_d[!(ls02_p_d$x==0),] 
ls05_p_d <- ls05_p_d[!(ls05_p_d$x==0),] 
ls07_p_d <- ls07_p_d[!(ls07_p_d$x==0),] 
ls08_p_d <- ls08_p_d[!(ls08_p_d$x==0),] 
ls09_p_d <- ls09_p_d[!(ls09_p_d$x==0),] 


#Scale y axis by multiplying by 10^4
ca01_p_d$y <- ca01_p_d$y*10^4 
ca02_p_d$y <- ca02_p_d$y*10^4 
ca03_p_d$y <- ca03_p_d$y*10^4 
ca04_p_d$y <- ca04_p_d$y*10^4 
ch01_p_d$y <- ch01_p_d$y*10^4 
ch02_p_d$y <- ch02_p_d$y*10^4 
ch03_p_d$y <- ch03_p_d$y*10^4 
ch04_p_d$y <- ch04_p_d$y*10^4 
ck01_p_d$y <- ck01_p_d$y*10^4 
ck02_p_d$y <- ck02_p_d$y*10^4 
ck03_p_d$y <- ck03_p_d$y*10^4 
ck04_p_d$y <- ck04_p_d$y*10^4 
cn01_p_d$y <- cn01_p_d$y*10^4 
cn02_p_d$y <- cn02_p_d$y*10^4 
ls01_p_d$y <- ls01_p_d$y*10^4 
ls02_p_d$y <- ls02_p_d$y*10^4 
ls05_p_d$y <- ls05_p_d$y*10^4 
ls07_p_d$y <- ls07_p_d$y*10^4 
ls08_p_d$y <- ls08_p_d$y*10^4 
ls09_p_d$y <- ls09_p_d$y*10^4 


########################################################
#                       SMC++                          #
########################################################

# set plotting values for SMC++
ca_s_file <- "./SMC++_MASKED/CA_7.26e-9_t10-10000.gen.csv"
ch_s_file <- "./SMC++_MASKED/CH_7.26e-9_t10-10000.gen.csv"
ck_s_file <- "./SMC++_MASKED/CK_7.26e-9_t10-10000.gen.csv"
cn_s_file <- "./SMC++_MASKED/CN_7.26e-9_t10-10000.gen.csv"
lean_s_file <- "./SMC++_MASKED/SN-lean_7.26e-9_t10-10000.gen.csv"
sis_s_file <- "./SMC++_MASKED/SN-siscowet_7.26e-9_t10-10000.gen.csv"


#read in data
ca_s_d <- read.csv(ca_s_file, header = TRUE)
ch_s_d <- read.csv(ch_s_file, header = TRUE)
ck_s_d<- read.csv(ck_s_file, header = TRUE)
cn_s_d <- read.csv(cn_s_file, header = TRUE)
lean_s_d <- read.csv(lean_s_file, header =TRUE)
sis_s_d <- read.csv(sis_s_file, header =TRUE)


#eliminate values with 0 time points
ca_s_d <- ca_s_d[!(ca_s_d$x==0),] 
ch_s_d <- ch_s_d[!(ch_s_d$x==0),] 
ck_s_d <- ck_s_d[!(ck_s_d$x==0),]
cn_s_d <- cn_s_d[!(cn_s_d$x==0),] 
lean_s_d <- lean_s_d[!(lean_s_d$x==0),] 
sis_s_d <- sis_s_d[!(sis_s_d$x==0),] 

# scale x axis based on generation time
gen_timeI <- 6

ca_s_d$x_year <- ca_s_d$x * gen_timeI
ch_s_d$x_year <- ch_s_d$x * gen_timeI
ck_s_d$x_year <- ck_s_d$x * gen_timeI
cn_s_d$x_year <- cn_s_d$x * gen_timeI

gen_timeII <- 16

lean_s_d$x_year <- lean_s_d$x * gen_timeII 
sis_s_d$x_year <- sis_s_d$x * gen_timeII






#colors <- c("#EA5869", "#FAD000", "#000080", "#000000", "#8b0000", "#FF0000")
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

a <- ggplot() +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = ck01_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ck02_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ck03_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ck04_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = cn01_p_d, aes(x=x, y=y, color=label), color = bsColor4, linewidth=1.25, alpha = 0.80) +
  geom_line(data = cn02_p_d, aes(x=x, y=y, color=label), color = bsColor4, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch01_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch02_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch03_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch04_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca01_p_d, aes(x=x, y=y, color=label), color =bsColor1, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca02_p_d, aes(x=x, y=y, color=label), color =bsColor1, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca03_p_d, aes(x=x, y=y, color=label), color =bsColor1, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca04_p_d, aes(x=x, y=y, color=label), color =bsColor1, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls01_p_d, aes(x=x, y=y, color=label), color =bsColor11, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls02_p_d, aes(x=x, y=y, color=label), color =bsColor11, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls07_p_d, aes(x=x, y=y, color=label), color =bsColor11, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls05_p_d, aes(x=x, y=y, color=label), color = bsColor14, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls08_p_d, aes(x=x, y=y, color=label), color = bsColor14, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ls09_p_d, aes(x=x, y=y, color=label), color = bsColor14, linewidth=1.25, alpha = 0.80) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  #scale_color_manual(name="Species", values=colors) +
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


b <- ggplot() +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = ck_s_d, aes(x=x_year, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = cn_s_d, aes(x=x_year, y=y, color=label), color = bsColor4, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch_s_d, aes(x=x_year, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca_s_d, aes(x=x_year, y=y, color=label), color =bsColor1, linewidth=1.25, alpha = 0.80) +
  geom_line(data = lean_s_d, aes(x=x_year, y=y, color=label), color =bsColor11, linewidth=1.25, alpha = 0.80) +
  geom_line(data = sis_s_d, aes(x=x_year, y=y, color=label), color =bsColor14, linewidth=1.25, alpha = 0.80) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  #scale_color_manual(name="Species", values=colors) +
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



plots <- plot_grid(a,b)
plots


wdII <- ("C:/Users/NJCB/Documents/Buffalo/Coregonus/Coregonus_PSMC-SMC_v2/figures_and_tables/SupFig7_Coregonus_Snamaycush_PSMC_SMC_no_bs/")
ggsave(paste(wdII, "./", "SupFig7_Coregonus_Snamaycush_PSMC_SMC_nobs_11-21-23", ".png", sep=""), plot=plots, units = "in", dpi=1000, height=4, width=8)
ggsave(paste(wdII, "./", "SupFig7_Coregonus_Snamaycush_PSMC_SMC_nobs_11-21-23", ".svg", sep=""), plot=plots, units = "in", dpi=1000, height=4, width=8)

