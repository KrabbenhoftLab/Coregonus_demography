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
ca01_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA01_bs.txt"
ca01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA01_main.txt"
ca02_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA02_bs.txt"
ca02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA02_main.txt"
ca03_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA03_bs.txt"
ca03_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA03_main.txt"
ca04_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA04_bs.txt"
ca04_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CA04_main.txt"
ch01_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH01_bs.txt"
ch01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH01_main.txt"
ch02_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH02_bs.txt"
ch02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH02_main.txt"
ch03_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH03_bs.txt"
ch03_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH03_main.txt"
ch04_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH04_bs.txt"
ch04_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CH04_main.txt"
ck01_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK01_bs.txt"
ck01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK01_main.txt"
ck02_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK02_bs.txt"
ck02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK02_main.txt"
ck03_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK03_bs.txt"
ck03_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK03_main.txt"
ck04_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK04_bs.txt"
ck04_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CK04_main.txt"
cn01_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CN01_bs.txt"
cn01_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CN01_main.txt"
cn02_p_bs_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CN02_bs.txt"
cn02_p_file <- "./PSMC_MASKED/PSMC_MASKED_param7_bootstrap_raw_bam/CN02_main.txt"


# read in PSMC data
ca01_p_bs <- read.delim(ca01_p_bs_file, header = TRUE, sep = "\t")
ca01_p_d <- read.delim(ca01_p_file, header = TRUE, sep = "\t")
ca02_p_bs <- read.delim(ca02_p_bs_file, header = TRUE, sep = "\t")
ca02_p_d <- read.delim(ca02_p_file, header = TRUE, sep = "\t")
ca03_p_bs <- read.delim(ca03_p_bs_file, header = TRUE, sep = "\t")
ca03_p_d <- read.delim(ca03_p_file, header = TRUE, sep = "\t")
ca04_p_bs <- read.delim(ca04_p_bs_file, header = TRUE, sep = "\t")
ca04_p_d <- read.delim(ca04_p_file, header = TRUE, sep = "\t")
ch01_p_bs <- read.delim(ch01_p_bs_file, header = TRUE, sep = "\t")
ch01_p_d <- read.delim(ch01_p_file, header = TRUE, sep = "\t")
ch02_p_bs <- read.delim(ch02_p_bs_file, header = TRUE, sep = "\t")
ch02_p_d <- read.delim(ch02_p_file, header = TRUE, sep = "\t")
ch03_p_bs <- read.delim(ch03_p_bs_file, header = TRUE, sep = "\t")
ch03_p_d <- read.delim(ch03_p_file, header = TRUE, sep = "\t")
ch04_p_bs <- read.delim(ch04_p_bs_file, header = TRUE, sep = "\t")
ch04_p_d <- read.delim(ch04_p_file, header = TRUE, sep = "\t")
ck01_p_bs <- read.delim(ck01_p_bs_file, header = TRUE, sep = "\t")
ck01_p_d <- read.delim(ck01_p_file, header = TRUE, sep = "\t")
ck02_p_bs <- read.delim(ck02_p_bs_file, header = TRUE, sep = "\t")
ck02_p_d <- read.delim(ck02_p_file, header = TRUE, sep = "\t")
ck03_p_bs <- read.delim(ck03_p_bs_file, header = TRUE, sep = "\t")
ck03_p_d <- read.delim(ck03_p_file, header = TRUE, sep = "\t")
ck04_p_bs <- read.delim(ck04_p_bs_file, header = TRUE, sep = "\t")
ck04_p_d <- read.delim(ck04_p_file, header = TRUE, sep = "\t")
cn01_p_bs <- read.delim(cn01_p_bs_file, header = TRUE, sep = "\t")
cn01_p_d <- read.delim(cn01_p_file, header = TRUE, sep = "\t")
cn02_p_bs <- read.delim(cn02_p_bs_file, header = TRUE, sep = "\t")
cn02_p_d <- read.delim(cn02_p_file, header = TRUE, sep = "\t")


# remove rows with "0" timepoints
ca01_p_bs <- ca01_p_bs[!(ca01_p_bs$x==0),] 
ca01_p_d <- ca01_p_d[!(ca01_p_d$x==0),] 
ca02_p_bs <- ca02_p_bs[!(ca02_p_bs$x==0),] 
ca02_p_d <- ca02_p_d[!(ca02_p_d$x==0),] 
ca03_p_bs <- ca03_p_bs[!(ca03_p_bs$x==0),] 
ca03_p_d <- ca03_p_d[!(ca03_p_d$x==0),] 
ca04_p_bs <- ca04_p_bs[!(ca04_p_bs$x==0),] 
ca04_p_d <- ca04_p_d[!(ca04_p_d$x==0),] 
ch01_p_bs <- ch01_p_bs[!(ch01_p_bs$x==0),] 
ch01_p_d <- ch01_p_d[!(ch01_p_d$x==0),] 
ch02_p_bs <- ch02_p_bs[!(ch02_p_bs$x==0),] 
ch02_p_d <- ch02_p_d[!(ch02_p_d$x==0),] 
ch03_p_bs <- ch03_p_bs[!(ch03_p_bs$x==0),] 
ch03_p_d <- ch03_p_d[!(ch03_p_d$x==0),] 
ch04_p_bs <- ch04_p_bs[!(ch04_p_bs$x==0),] 
ch04_p_d <- ch04_p_d[!(ch04_p_d$x==0),] 
ck01_p_bs <- ck01_p_bs[!(ck01_p_bs$x==0),] 
ck01_p_d <- ck01_p_d[!(ck01_p_d$x==0),] 
ck02_p_bs <- ck02_p_bs[!(ck02_p_bs$x==0),] 
ck02_p_d <- ck02_p_d[!(ck02_p_d$x==0),] 
ck03_p_bs <- ck03_p_bs[!(ck03_p_bs$x==0),] 
ck03_p_d <- ck03_p_d[!(ck03_p_d$x==0),] 
ck04_p_bs <- ck04_p_bs[!(ck04_p_bs$x==0),] 
ck04_p_d <- ck04_p_d[!(ck04_p_d$x==0),] 
cn01_p_bs <- cn01_p_bs[!(cn01_p_bs$x==0),] 
cn01_p_d <- cn01_p_d[!(cn01_p_d$x==0),] 
cn02_p_bs <- cn02_p_bs[!(cn02_p_bs$x==0),] 
cn02_p_d <- cn02_p_d[!(cn02_p_d$x==0),] 

#Scale y axis by multiplying by 10^4
ca01_p_bs$y <- ca01_p_bs$y*10^4 
ca01_p_d$y <- ca01_p_d$y*10^4 
ca02_p_bs$y <- ca02_p_bs$y*10^4 
ca02_p_d$y <- ca02_p_d$y*10^4 
ca03_p_bs$y <- ca03_p_bs$y*10^4 
ca03_p_d$y <- ca03_p_d$y*10^4 
ca04_p_bs$y <- ca04_p_bs$y*10^4 
ca04_p_d$y <- ca04_p_d$y*10^4 
ch01_p_bs$y <- ch01_p_bs$y*10^4 
ch01_p_d$y <- ch01_p_d$y*10^4 
ch02_p_bs$y <- ch02_p_bs$y*10^4 
ch02_p_d$y <- ch02_p_d$y*10^4 
ch03_p_bs$y <- ch03_p_bs$y*10^4 
ch03_p_d$y <- ch03_p_d$y*10^4 
ch04_p_bs$y <- ch04_p_bs$y*10^4 
ch04_p_d$y <- ch04_p_d$y*10^4 
ck01_p_bs$y <- ck01_p_bs$y*10^4 
ck01_p_d$y <- ck01_p_d$y*10^4 
ck02_p_bs$y <- ck02_p_bs$y*10^4 
ck02_p_d$y <- ck02_p_d$y*10^4 
ck03_p_bs$y <- ck03_p_bs$y*10^4 
ck03_p_d$y <- ck03_p_d$y*10^4 
ck04_p_bs$y <- ck04_p_bs$y*10^4 
ck04_p_d$y <- ck04_p_d$y*10^4 
cn01_p_bs$y <- cn01_p_bs$y*10^4 
cn01_p_d$y <- cn01_p_d$y*10^4 
cn02_p_bs$y <- cn02_p_bs$y*10^4 
cn02_p_d$y <- cn02_p_d$y*10^4 

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
bsColor1 <- brightness("#EA5869", 0.9) # make brighter version of the base color (good only for black)
bsColor2 <- brightness("#FAD000", 0.9)
bsColor3 <- brightness("#000080", 0.9)
bsColor4 <- brightness("#000000", 0.1)
bsColor5 <- brightness("#B5B35C", 0.9)
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



#plot PSMC

a <- ggplot(ca01_p_d, aes(x=x, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = ca01_p_d, aes(x=x, y=y, color=label), linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca01_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor1, alpha=0.1)) +
  geom_line(data = ca02_p_d, aes(x=x, y=y, color=label), linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca02_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor1, alpha=0.1)) +
  geom_line(data = ca03_p_d, aes(x=x, y=y, color=label), linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca03_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor1, alpha=0.1)) +
  geom_line(data = ca04_p_d, aes(x=x, y=y, color=label), linewidth=1.25, alpha = 0.80) +
  geom_line(data = ca04_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor1, alpha=0.1)) +
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

c <- ggplot(ch01_p_d, aes(x=x, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = ch01_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch01_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor2, alpha=0.1)) +
  geom_line(data = ch02_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch02_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor2, alpha=0.1)) +
  geom_line(data = ch03_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch03_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor2, alpha=0.1)) +
  geom_line(data = ch04_p_d, aes(x=x, y=y, color=label), color = bsColor2, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ch04_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor2, alpha=0.1)) +
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

e <- ggplot(ck01_p_d, aes(x=x, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = ck01_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ck01_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor3, alpha=0.1)) +
  geom_line(data = ck02_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ck02_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor3, alpha=0.1)) +
  geom_line(data = ck03_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ck03_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor3, alpha=0.1)) +
  geom_line(data = ck04_p_d, aes(x=x, y=y, color=label), color = bsColor3, linewidth=1.25, alpha = 0.80) +
  geom_line(data = ck04_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor3, alpha=0.1)) +
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

e

g <- ggplot(cn01_p_d, aes(x=x, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  #annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
  #fill = "black", alpha=0.6) +
  geom_line(data = cn01_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor4, alpha=0.1)) +
  geom_line(data = cn02_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor4, alpha=0.1)) +
  geom_line(data = cn01_p_d, aes(x=x, y=y, color=label), color = c("#000000"), linewidth=1.25, alpha = 0.80) +
  geom_line(data = cn02_p_d, aes(x=x, y=y, color=label), color = c("#000000"), linewidth=1.25, alpha = 0.80) +
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

g


#plot SMC++

b <- ggplot(ca_s_d, aes(x=x_year, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(linewidth = 1.25, color = bsColor1) +
  geom_line(data = ca_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor1, alpha=0.2)) +
  geom_line(data = ca_s_d, aes(x=x_year, y=y, group=label), color=bsColor1) +
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

d <- ggplot(ch_s_d, aes(x=x_year, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(linewidth = 1.25, color = bsColor2) +
  geom_line(data = ch_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor2, alpha=0.2)) +
  geom_line(data = ch_s_d, aes(x=x_year, y=y, group=label), color=bsColor2) +
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
  
f <- ggplot(ck_s_d, aes(x=x_year, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(linewidth = 1.25, color = bsColor3) +
  geom_line(data = ck_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor3, alpha=0.1)) +
  geom_line(data = ck_s_d, aes(x=x_year, y=y, group=label), color=bsColor3) +
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
  
f

h <- ggplot(cn_s_d, aes(x=x_year, y=y, color=label)) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(linewidth = 1.25, color = bsColor4) +
  geom_line(data = cn_s_bs, aes(x=x_year, y=y, group=label), color=alpha(bsColor4, alpha=0.1)) +
  geom_line(data = cn_s_d, aes(x=x_year, y=y, group=label), color=bsColor4) +
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

h

plots <- plot_grid(a,
                   b,
                   c,
                   d,
                   e,
                   f,
                   g,
                   h,
                   ncol = 2)

  
plots  
  

wdII <- ("C:/Users/NJCB/Documents/Buffalo/Coregonus/Coregonus_PSMC-SMC_v2/figures_and_tables/SupFig2-Coregonus_PSMC_SMC_bs_plot_grid/")
ggsave(paste(wdII, "./", "SupFig2_Coregonus_PSMC_SMC_bs_plot_grid_11-20-23", ".png", sep=""), plot=plots, units = "in", dpi=1000, height=10, width=8)
ggsave(paste(wdII, "./", "SupFig2_Coregonus_PSMC_SMC_bs_plot_grid_11-20-23", ".svg", sep=""), plot=plots, units = "in", dpi=1000, height=10, width=8)
