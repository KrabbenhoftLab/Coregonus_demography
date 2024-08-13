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
#                       SMC++                          #
########################################################

# set plotting values for SMC++
ca_s_file <- "./SMC++_MASKED/CA_7.26e-9_t10-10000.gen.csv"
ca_s_bs_file <- "./SMC++_MASKED/CA_bs_all.gen.csv"
ch_s_file <- "./SMC++_MASKED/CH_7.26e-9_t10-10000.gen.csv"
ch_s_bs_file <- "./SMC++_MASKED/CH_bs_all.gen.csv"
ck_s_file <- "./SMC++_MASKED/CK_7.26e-9_t10-10000.gen.csv"
ck_s_bs_file <- "./SMC++_MASKED/CK_bs_all.gen.csv"
cn_s_file <- "./SMC++_MASKED/CN_7.26e-9_t10-10000.gen.csv"
cn_s_bs_file <- "./SMC++_MASKED/CN_bs_all.gen.csv"


#read in data
ca_s_d <- read.csv(ca_s_file, header = TRUE)
ca_s_bs <- read.csv(ca_s_bs_file, header = TRUE)
ch_s_d <- read.csv(ch_s_file, header = TRUE)
ch_s_bs <- read.csv(ch_s_bs_file, header = TRUE)
ck_s_d<- read.csv(ck_s_file, header = TRUE)
ck_s_bs <- read.csv(ck_s_bs_file, header = TRUE)
cn_s_d <- read.csv(cn_s_file, header = TRUE)
cn_s_bs <- read.csv(cn_s_bs_file, header = TRUE) 

#eliminate values with 0 time points
ca_s_d <- ca_s_d[!(ca_s_d$x==0),] 
ca_s_bs <- ca_s_bs[!(ca_s_bs$x==0),] 
ch_s_d <- ch_s_d[!(ch_s_d$x==0),] 
ch_s_bs <- ch_s_bs[!(ch_s_bs$x==0),] 
ck_s_d <- ck_s_d[!(ck_s_d$x==0),]
ck_s_bs <- ck_s_bs[!(ck_s_bs$x==0),] 
cn_s_d <- cn_s_d[!(cn_s_d$x==0),] 
cn_s_bs <- cn_s_bs[!(cn_s_bs$x==0),] 

# scale x axis based on generation time
gen_time <- 6

ca_s_d$x_year <- ca_s_d$x * gen_time
ca_s_bs$x_year <- ca_s_bs$x * gen_time
ch_s_d$x_year <- ch_s_d$x * gen_time
ch_s_bs$x_year <- ch_s_bs$x * gen_time
ck_s_d$x_year <- ck_s_d$x * gen_time
ck_s_bs$x_year <- ck_s_bs$x * gen_time
cn_s_d$x_year <- cn_s_d$x * gen_time
cn_s_bs$x_year <- cn_s_bs$x * gen_time


# color blind safe palette
# color order matches alphabetical order of species
# i.e. CA, CH, CK, CN, CZ
colors <- c("#EA5869", "#EA5869", "#EA5869", "#EA5869", "#FAD000", "#FAD000", "#FAD000", "#FAD000", "#000080", "#000080", "#000080", "#000080", "#000000", "#000000", "#8b0000", "#8b0000", "#FF0000", "#8b0000", "#FF0000", "#FF0000")
bsColor1 <- c("#EA5869") # make brighter version of the base color (good only for black)
bsColor2 <- c("#FAD000")
bsColor3 <- c("#000080")
bsColor4 <- c("#000000")
bsColor5 <- c("#B5B35C")
bsColor6 <- c("#2ee009") #*
bsColor7 <- c("#29d209")
bsColor8 <- c("#24c509")
bsColor9 <- c("#20b808")
bsColor10 <- c("#1bab08")
bsColor11 <- c("#139107") #*
bsColor12 <- c("#0c7805")
bsColor13 <- c("#056003")
bsColor14 <- c("#024902") #*

# geologic time boundaries, from ICS v2020/03 (https://stratigraphy.org/chart)
pliestocene <- c(11700, Inf) # use Inf to plot properly
# last glacial maximum, from Clark et al. 2009, https://science.sciencemag.org/content/325/5941/710.abstract
LGM <- c(19000, 33000)
# ice retreat from Great Lakes, from Dalton et al. 2020 https://www.sciencedirect.com/science/article/pii/S0277379119307619
deglaciation <- c(13500, 15500)


#plot SMC++

sf3 <- ggplot(ca_s_d, aes(x=x_year, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(linewidth = 1.25, color = bsColor1) +
  geom_line(data = cn_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor4, alpha=0.1)) +
  geom_line(data = cn_s_d, aes(x=x_year, y=y, group=label), color=bsColor4, linewidth = 1.25) +
  geom_line(data = ck_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor3, alpha=0.1)) +
  geom_line(data = ck_s_d, aes(x=x_year, y=y, group=label), color=bsColor3, linewidth = 1.25) +
  geom_line(data = ca_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor1, alpha=0.1)) +
  geom_line(data = ca_s_d, aes(x=x_year, y=y, group=label), color=bsColor1, linewidth = 1.25) +
  geom_line(data = ch_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor2, alpha=0.1)) +
  geom_line(data = ch_s_d, aes(x=x_year, y=y, group=label), color=bsColor2, linewidth = 1.25) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL#,
    #limits=c(1000, 2800000)
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

sf3


wdII <- ("C:/Users/NJCB/Documents/Buffalo/Coregonus/Coregonus_PSMC-SMC_v2/figures_and_tables/SupFig3-Coregonus_SMC_bootstraps/")
ggsave(paste(wdII, "./", "SupFig3_Coregonus_SMC_bs_11-20-23", ".png", sep=""), plot=sf3, units = "in", dpi=1000, height=6, width=8)
ggsave(paste(wdII, "./", "SupFig3_Coregonus_SMC_bs_11-20-23", ".svg", sep=""), plot=sf3, units = "in", dpi=1000, height=6, width=8)







