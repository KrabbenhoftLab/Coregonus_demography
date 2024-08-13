#################################################################
# Script to perform PCA analyses with genome-wide SNP data
# Dan MacGuigan
#################################################################

library(caret)
library(ggplot2)
library(ggfortify)
library(ape)
library(phytools)
library(adegenet)
library(reshape2)
library(gridBase)
library(grid)
library(gridExtra)
library(Cairo)
library(cowplot)
library(vcfR)
library(plyr)
library(ggrepel)
library(svglite)

# set working directory
setwd("H:/KrabbLab/Cisco_demography/adegenet_revisions/")

# read in VCF
# may take a few minutes
vcf <- read.vcfR("./data/Coregonus_PCA_biallelic_no_missing.maf0.05.vcf.gz")

# convert to genlight format
# may take a few minutes
d <- vcfR2genlight(vcf)

# read in sample data
# one row per sample, sample info (species, ID, location, etc) in columns
# sample (row) order needs to match sample order in VCF file
data <- read.csv("./data/sampleInfo.csv")

# perform PCA
# may take a few minutes
pca <- glPca(d, nf = 5)

# save pca object to file for later replotting
saveRDS(pca, "./data/pca.maf0.05.rds")

# reload pca object if needed
#pca <- readRDS("./data/pca.maf0.05.rds")

# quick scatter plot of PC axis 1 vs 2
scatter(pca)

# function to plot PC axes
plot_genlight_PCA <- function(df, group_col, sample_col, show_sample_names, show_legend, group_names, group_colors, x.axis, y.axis){
  find_hull <- function(df) df[chull(df[,x.axis], df[,y.axis]), ] # function to get convex hull
  hulls <- ddply(df, group_col, find_hull) 
  
  if(show_sample_names && show_legend){
    p <- ggplot(df, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""), colour=group_col)) +
      geom_polygon(data=hulls, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""),
                                          color=group_col, fill=group_col), alpha=0.2) +
      scale_colour_manual(values = group_colors, labels = group_names, name = "Species") +
      geom_point(cex=5, alpha=0.7) +
      scale_fill_manual(values = group_colors, labels = group_names, name = "Species") +
      theme_minimal() +
      theme(legend.position="right", legend.box = "horizontal") +
      xlab(paste("PC ", x.axis, " (", round(100*pca$eig[x.axis]/sum(pca$eig), digits = 2), "% variation)", sep = "")) +
      ylab(paste("PC ", y.axis, " (", round(100*pca$eig[y.axis]/sum(pca$eig), digits = 2), "% variation)", sep = "")) +
      geom_label_repel(aes_string(label = sample_col),
                       box.padding   = 0.35, 
                       point.padding = 0.5,
                       segment.color = 'grey50', show.legend = FALSE) 
  } else if(show_sample_names && !show_legend) {
    p <- ggplot(df, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""), colour=group_col)) +
      geom_polygon(data=hulls, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""),
                                          color=group_col, fill=group_col), alpha=0.2) +
      scale_colour_manual(values = group_colors, labels = group_names, name = "Species") +
      geom_point(cex=5, alpha=0.7) +
      scale_fill_manual(values = group_colors, labels = group_names, name = "Species") +
      theme_minimal() +
      theme(legend.position="none", legend.box = "horizontal") +
      xlab(paste("PC ", x.axis, " (", round(100*pca$eig[x.axis]/sum(pca$eig), digits = 2), "% variation)", sep = "")) +
      ylab(paste("PC ", y.axis, " (", round(100*pca$eig[y.axis]/sum(pca$eig), digits = 2), "% variation)", sep = "")) +
      geom_label_repel(aes_string(label = sample_col),
                       box.padding   = 0.35, 
                       point.padding = 0.5,
                       segment.color = 'grey50', show.legend = FALSE) 
    } else if(show_legend) {
    p <- ggplot(df, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""), colour=group_col)) +
      geom_polygon(data=hulls, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""),
                                          color=group_col, fill=group_col), alpha=0.2) +
      scale_colour_manual(values = group_colors, labels = group_names, name = "Species") +
      geom_point(cex=5, alpha=0.7) +
      scale_fill_manual(values = group_colors, labels = group_names, name = "Species") +
      theme_minimal() +
      theme(legend.position="right", legend.box = "horizontal") +
      xlab(paste("PC ", x.axis, " (", round(100*pca$eig[x.axis]/sum(pca$eig), digits = 2), "% variation)", sep = "")) +
      ylab(paste("PC ", y.axis, " (", round(100*pca$eig[y.axis]/sum(pca$eig), digits = 2), "% variation)", sep = ""))
  } else {
    p <- ggplot(df, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""), colour=group_col)) +
      geom_polygon(data=hulls, aes_string(x=paste("PC", x.axis, sep=""), y=paste("PC", y.axis, sep=""),
                                          color=group_col, fill=group_col), alpha=0.2) +
      scale_colour_manual(values = group_colors, labels = group_names, name = "Species") +
      geom_point(cex=5, alpha=0.7) +
      scale_fill_manual(values = group_colors, labels = group_names, name = "Species") +
      theme_minimal() +
      theme(legend.position="none", legend.box = "horizontal") +
      xlab(paste("PC ", x.axis, " (", round(100*pca$eig[x.axis]/sum(pca$eig), digits = 2), "% variation)", sep = "")) +
      ylab(paste("PC ", y.axis, " (", round(100*pca$eig[y.axis]/sum(pca$eig), digits = 2), "% variation)", sep = ""))
  }
  return(p)
}

# make plotting data frame
points <- data.frame(pca$scores)
points$species <- data$species # add some other info like species
points$sample <- data$sample # sample ID
points <- points[order(points$species), ] # order data by species

# vector of species names and colors
species <- c(bquote(italic("C. artedi")),
             bquote(italic("C. hoyi")), 
             bquote(italic("C. kiyi")), 
             bquote(italic("C. nigripinnis"))) 
colors <- c("#EA5869", "#fad000", "#000080", "#000000")

# quick test plot to make sure everything is working
test <- plot_genlight_PCA(df = points,
                          group_col = "species",
                          sample_col = "sample",
                          show_sample_names = TRUE,
                          show_legend = TRUE,
                          group_names = species,
                          group_colors = colors,
                          x.axis = 1,
                          y.axis = 2)
test



# make a bunch of ggplot objects comparing different PC axes
plot_pc1_pc2 <- plot_genlight_PCA(df = points,
                                  group_col = "species",
                                  sample_col = "sample",
                                  show_sample_names = TRUE,
                                  show_legend = FALSE,
                                  group_names = species,
                                  group_colors = colors,
                                  x.axis = 1,
                                  y.axis = 2)

plot_pc1_pc3 <- plot_genlight_PCA(df = points,
                                  group_col = "species",
                                  sample_col = "sample",
                                  show_sample_names = TRUE,
                                  show_legend = FALSE,
                                  group_names = species,
                                  group_colors = colors,
                                  x.axis = 1,
                                  y.axis = 3)

plot_pc1_pc4 <- plot_genlight_PCA(df = points,
                                  group_col = "species",
                                  sample_col = "sample",
                                  show_sample_names = TRUE,
                                  show_legend = FALSE,
                                  group_names = species,
                                  group_colors = colors,
                                  x.axis = 1,
                                  y.axis = 4)

plot_pc2_pc3 <- plot_genlight_PCA(df = points,
                                  group_col = "species",
                                  sample_col = "sample",
                                  show_sample_names = TRUE,
                                  show_legend = FALSE,
                                  group_names = species,
                                  group_colors = colors,
                                  x.axis = 2,
                                  y.axis = 3)

plot_pc2_pc4 <- plot_genlight_PCA(df = points,
                                  group_col = "species",
                                  sample_col = "sample",
                                  show_sample_names = TRUE,
                                  show_legend = FALSE,
                                  group_names = species,
                                  group_colors = colors,
                                  x.axis = 2,
                                  y.axis = 4)

plot_pc3_pc4 <- plot_genlight_PCA(df = points,
                                  group_col = "species",
                                  sample_col = "sample",
                                  show_sample_names = TRUE,
                                  show_legend = FALSE,
                                  group_names = species,
                                  group_colors = colors,
                                  x.axis = 3,
                                  y.axis = 4)

temp  <- plot_genlight_PCA(df = points,
                           group_col = "species",
                           sample_col = "sample",
                           show_sample_names = TRUE,
                           show_legend = TRUE,
                           group_names = species,
                           group_colors = colors,
                           x.axis = 3,
                           y.axis = 4)
l <- get_legend(temp) # extract legend because we only want to plot it once
ggdraw(l) # check to make sure the legend looks ok

# grid plot
pdf("./PCA/figures_maf0.05/PCA.pdf", height=10,width=15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,3))) # specify blank 8x3 grid with set widths
vplayout <- function(x,y) viewport(layout.pos.row = x, layout.pos.col = y)
print(ggdraw(l), vp = vplayout(1,3)) # plot legend
print(plot_pc1_pc2, vp = vplayout(1,1))
print(plot_pc1_pc3, vp = vplayout(2,1))
print(plot_pc1_pc4, vp = vplayout(3,1))
print(plot_pc2_pc3, vp = vplayout(2,2))
print(plot_pc2_pc4, vp = vplayout(3,2))
print(plot_pc3_pc4, vp = vplayout(3,3))
dev.off()

svglite("./PCA/figures_maf0.05/PCA.svg", height=10,width=15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,3))) # specify blank 8x3 grid with set widths
vplayout <- function(x,y) viewport(layout.pos.row = x, layout.pos.col = y)
print(ggdraw(l), vp = vplayout(1,3)) # plot legend
print(plot_pc1_pc2, vp = vplayout(1,1))
print(plot_pc1_pc3, vp = vplayout(2,1))
print(plot_pc1_pc4, vp = vplayout(3,1))
print(plot_pc2_pc3, vp = vplayout(2,2))
print(plot_pc2_pc4, vp = vplayout(3,2))
print(plot_pc3_pc4, vp = vplayout(3,3))
dev.off()

# plot just PC1 vs PC2
pdf("./PCA/figures_maf0.05/PCA_PC1-2.pdf", height=5, width=8)
plot_genlight_PCA(df = points,
                  group_col = "species",
                  sample_col = "sample",
                  show_sample_names = TRUE,
                  show_legend = TRUE,
                  group_names = species,
                  group_colors = colors,
                  x.axis = 1,
                  y.axis = 2)
dev.off()

svglite("./PCA/figures_maf0.05/PCA_PC1-2.svg", height=5, width=8)
plot_genlight_PCA(df = points,
                  group_col = "species",
                  sample_col = "sample",
                  show_sample_names = TRUE,
                  show_legend = TRUE,
                  group_names = species,
                  group_colors = colors,
                  x.axis = 1,
                  y.axis = 2)
dev.off()

