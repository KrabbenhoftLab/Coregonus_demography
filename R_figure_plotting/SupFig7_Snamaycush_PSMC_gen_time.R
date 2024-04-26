##############################################################################################################################

# script to plot PSMC & SMC++ Snamaycush multi-generation time
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
library(extrafont)

# specify the following parameters
wd <- "C:/Users/NJCB/Documents/Buffalo/Coregonus/Coregonus_PSMC-SMC_v2/results/Snamaycush_PSMC_MASKED_gen_time/" # top level working directory
setwd(wd)
getwd()
list.files()


########################################################
#                       PSMC                           #
########################################################
# set plotting values for PSMC
#These samples are all LS01

g6_p_bs_file <- "./LS01_g6_bs.txt"
g6_p_file <- "./LS01_g6_main.txt"
g8_p_bs_file <- "./LS01_g8_bs.txt"
g8_p_file <- "./LS01_g8_main.txt"
g10_p_bs_file <- "./LS01_g10_bs.txt"
g10_p_file <- "./LS01_g10_main.txt"
g12_p_bs_file <- "./LS01_g12_bs.txt"
g12_p_file <- "./LS01_g12_main.txt"
g14_p_bs_file <- "./LS01_g14_bs.txt"
g14_p_file <- "./LS01_g14_main.txt"
g16_p_bs_file <- "./LS01_g16_bs.txt"
g16_p_file <- "./LS01_g16_main.txt"
g18_p_bs_file <- "./LS01_g18_bs.txt"
g18_p_file <- "./LS01_g18_main.txt"
g20_p_bs_file <- "./LS01_g20_bs.txt"
g20_p_file <- "./LS01_g20_main.txt"



# read in PSMC data
g6_p_bs <- read.delim(g6_p_bs_file, header = TRUE, sep = "\t") 
g6_p_d <- read.delim(g6_p_file, header = TRUE, sep = "\t") 
g8_p_bs <- read.delim(g8_p_bs_file, header = TRUE, sep = "\t") 
g8_p_d <- read.delim(g8_p_file, header = TRUE, sep = "\t") 
g10_p_bs <- read.delim(g10_p_bs_file, header = TRUE, sep = "\t") 
g10_p_d <- read.delim(g10_p_file, header = TRUE, sep = "\t") 
g12_p_bs <- read.delim(g12_p_bs_file, header = TRUE, sep = "\t") 
g12_p_d <- read.delim(g12_p_file, header = TRUE, sep = "\t") 
g14_p_bs <- read.delim(g14_p_bs_file, header = TRUE, sep = "\t") 
g14_p_d <- read.delim(g14_p_file, header = TRUE, sep = "\t") 
g16_p_bs <- read.delim(g16_p_bs_file, header = TRUE, sep = "\t") 
g16_p_d <- read.delim(g16_p_file, header = TRUE, sep = "\t") 
g18_p_bs <- read.delim(g18_p_bs_file, header = TRUE, sep = "\t") 
g18_p_d <- read.delim(g18_p_file, header = TRUE, sep = "\t") 
g20_p_bs <- read.delim(g20_p_bs_file, header = TRUE, sep = "\t") 
g20_p_d <- read.delim(g20_p_file, header = TRUE, sep = "\t") 


# remove rows with "0" timepoints
g6_p_bs <- g6_p_bs[!(g6_p_bs$x==0),]
g6_p_d <- g6_p_d[!(g6_p_d$x==0),]
g8_p_bs <- g8_p_bs[!(g8_p_bs$x==0),]
g8_p_d <- g8_p_d[!(g8_p_d$x==0),]
g10_p_bs <- g10_p_bs[!(g10_p_bs$x==0),]
g10_p_d <- g10_p_d[!(g10_p_d$x==0),]
g12_p_bs <- g12_p_bs[!(g12_p_bs$x==0),]
g12_p_d <- g12_p_d[!(g12_p_d$x==0),]
g14_p_bs <- g14_p_bs[!(g14_p_bs$x==0),]
g14_p_d <- g14_p_d[!(g14_p_d$x==0),]
g16_p_bs <- g16_p_bs[!(g16_p_bs$x==0),]
g16_p_d <- g16_p_d[!(g16_p_d$x==0),]
g18_p_bs <- g18_p_bs[!(g18_p_bs$x==0),]
g18_p_d <- g18_p_d[!(g18_p_d$x==0),]
g20_p_bs <- g20_p_bs[!(g20_p_bs$x==0),]
g20_p_d <- g20_p_d[!(g20_p_d$x==0),]

#Scale y axis by multiplying by 10^4
g6_p_bs$y <- g6_p_bs$y*10^4
g6_p_d$y <- g6_p_d$y*10^4
g8_p_bs$y <- g8_p_bs$y*10^4
g8_p_d$y <- g8_p_d$y*10^4
g10_p_bs$y <- g10_p_bs$y*10^4
g10_p_d$y <- g10_p_d$y*10^4
g12_p_bs$y <- g12_p_bs$y*10^4
g12_p_d$y <- g12_p_d$y*10^4
g14_p_bs$y <- g14_p_bs$y*10^4
g14_p_d$y <- g14_p_d$y*10^4
g16_p_bs$y <- g16_p_bs$y*10^4
g16_p_d$y <- g16_p_d$y*10^4
g18_p_bs$y <- g18_p_bs$y*10^4
g18_p_d$y <- g18_p_d$y*10^4
g20_p_bs$y <- g20_p_bs$y*10^4
g20_p_d$y <- g20_p_d$y*10^4



# color blind safe palette
# color order matches alphabetical order of species
# i.e. CA, CH, CK, CN, CZ
colors <- c("#8b0000", "#8b0000", "#8b0000")
colorsII <- c("#ff2929", "#e22026", "#c51722", "#a90f1e", "#8e0819", "#730413", "#59010c", "#410000")
bsColor1 <- c("#EA5869") # make brighter version of the base color (good only for black)
bsColor2 <- c("#FAD000")
bsColor3 <- c("#000080")
bsColor4 <- c("#000000")
bsColor6<- c("#ff2929")
bsColor8<- c("#e22026")
bsColor10<- c("#c51722")
bsColor12<- c("#a90f1e")
bsColor14<- c("#8e0819")
bsColor16<- c("#730413")
bsColor18<- c("#59010c")
bsColor20<- c("#410000")
# geologic time boundaries, from ICS v2020/03 (https://stratigraphy.org/chart)
pliestocene <- c(11700, Inf) # use Inf to plot properly
# last glacial maximum, from Clark et al. 2009, https://science.sciencemag.org/content/325/5941/710.abstract
LGM <- c(19000, 33000)
# ice retreat from Great Lakes, from Dalton et al. 2020 https://www.sciencedirect.com/science/article/pii/S0277379119307619
deglaciation <- c(13500, 15500)



#plot PSMC


a <- ggplot(g6_p_d, aes(x=x, y=y, color=bsColor6)) +
  #annotate("rect", xmin = pliestocene[1], xmax = pliestocene[2], ymin = 0, ymax = Inf,
  #         fill = "gray", alpha=0.6) +
  annotate("rect", xmin = LGM[1], xmax = LGM[2], ymin = 0, ymax = Inf,
           fill = "gray", alpha=0.6) +
  annotate("rect", xmin = deglaciation[1], xmax = deglaciation[2], ymin = 0, ymax = Inf,
           fill = "black", alpha=0.6) +
  geom_line(data = g6_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor6, alpha=0.15)) +
  geom_line(data = g8_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor8, alpha=0.15)) +
  geom_line(data = g10_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor10, alpha=0.15)) +
  geom_line(data = g12_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor12, alpha=0.15)) +
  geom_line(data = g14_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor14, alpha=0.15)) +
  geom_line(data = g16_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor16, alpha=0.15)) +
  geom_line(data = g18_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor18, alpha=0.15)) +
  geom_line(data = g20_p_bs, aes(x=x, y=y, group=label), color=alpha(bsColor20, alpha=0.15)) +
  geom_line(data = g6_p_d, aes(x=x, y=y), color = bsColor6, linewidth=1.25) +
  geom_line(data = g8_p_d, aes(x=x, y=y), color = bsColor8, linewidth=1.25) +
  geom_line(data = g10_p_d, aes(x=x, y=y), color = bsColor10, linewidth=1.25) +
  geom_line(data = g12_p_d, aes(x=x, y=y), color = bsColor12, linewidth=1.25) +
  geom_line(data = g14_p_d, aes(x=x, y=y), color = bsColor14, linewidth=1.25) +
  geom_line(data = g16_p_d, aes(x=x, y=y), color = bsColor16, linewidth=1.25) +
  geom_line(data = g18_p_d, aes(x=x, y=y), color = bsColor18, linewidth=1.25) +
  geom_line(data = g20_p_d, aes(x=x, y=y), color = bsColor20, linewidth=1.25) +
  ylab(bquote("")) +
  xlab(bquote("")) +
  scale_color_manual(name = "gen time", values = colorsII) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL,
    limits=c(1000, 3000000)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    minor_breaks = NULL
  ) +
  annotation_logticks(sides="bl") + # log scale tick marks for both axes
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  coord_cartesian(ylim=c(3000,1000000))


a


wdII <- ("C:/Users/NJCB/Documents/Buffalo/Coregonus/Coregonus_PSMC-SMC_v2/figures_and_tables/SupFig7_Snamaycush_PSMC_gen_time/")
ggsave(paste(wdII, "./", "SupFig7_Snamaycush_PSMC_gen_time_12-06-23", ".png", sep=""), plot=a, units = "in", dpi=1000, height=8, width=10)
ggsave(paste(wdII, "./", "SupFig7_Snamaycush_PSMC_gen_time_12-06-23", ".svg", sep=""), plot=a, units = "in", dpi=1000, height=8, width=10)
