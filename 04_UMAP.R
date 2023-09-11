########################################
#' Generate PBMC UMAP to visualize markers distribution and cell populations
#' First steps in analysis to have an overview of the data and make sure that everything looks normal
#' Karine Thai - 2023/06/08, Updated on 2023/07/24
########################################

######################
#' 0 - load required packages
######################

library(tidyverse)
library(ggplot2)
library(tools)
library(umap)
library(RColorBrewer)
library(tidyr)
library(ggplot2)

######################
#' 1 - Load list of tables of expression, normalized and clean from script 03
######################
load('Output/RData/ListTablesExpCut.RData')

######################
#' 2 - Check the min/max number of cells to make sure that there are enough cells for the downsample
######################

source("~/Documents/PhD/Coding/R_Functions/Min_number_of_cells_in_samples.R")
source("~/Documents/PhD/Coding/R_Functions/Max_number_of_cells_in_samples.R")
Min_number_of_cells_in_samples(ListTablesExpCut)
Max_number_of_cells_in_samples(ListTablesExpCut)

######################
#' 3 - Create a table containing a downsample of each sample
#' Will be used to generate UMAP
######################

source("~/Documents/PhD/Coding/R_Functions/Create_Table_UMAP_RandomEvents.R")

set.seed(1)
Table4UMAP_PBMC = TableExprsDownsample(ListTablesExpCut,1500)

tail(row.names(Table4UMAP_PBMC))
colnames(Table4UMAP_PBMC)


######################
#' 4 - Scale the data
######################

Scaled_Table4UMAP_PBMC=scale(Table4UMAP_PBMC[,1:24])

######################
#' 5 - Generate UMAP with scaled parameters
######################

Scaled_UMAP = umap(Scaled_Table4UMAP_PBMC[,c(1,3,5:11,13:24)],n_components=2,  n_neighbors=15,  min_dist=0.6)
#save(Scaled_UMAP, file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/RData/20230609_Scaled_UMAP.RData")

######################
#' 6 - Add UMAP coordinates to table
######################

df_UMAP_PBMC = data.frame(cbind(Table4UMAP_PBMC, Scaled_UMAP$layout))
colnames(df_UMAP_PBMC)[ncol(Table4UMAP_PBMC)+1:2]=c("UMAP1","UMAP2")
#save(df_UMAP_PBMC, file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/RData/20230609_df_UMAP_PBMC.RData")

######################
#' 7 - Plot and save UMAP
######################
pdf(file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/Plots/20230609_Scaled_UMAP.pdf",
    width = 12, 
    height = 9) 
ggplot(df_UMAP_PBMC,aes(UMAP1,UMAP2))+geom_point(size=0.2)+theme_bw()+theme(text = element_text(size=20))
dev.off()

######################
#' 8 - Plot and save markers on UMAP
######################

molecules=colnames(df_UMAP_PBMC)[c(1,3,5:24)]

plots=list()
for(mol in molecules){
  plots[[mol]]=ggplot(df_UMAP_PBMC,aes_string("UMAP1","UMAP2",col=mol))+geom_point(size=0.25)+theme_minimal(base_size = 8.5)+theme(text=element_text(size=14))+ scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(9,"RdYlBu")))
}

markers_on_umap = cowplot::plot_grid(plotlist=plots[c(2,3,16,9,4,7,15,14,11,17,8,13,6,18,5,19,20,21,12,22)], nrow=4)

pdf(file = "Output/Plots/20230609_MarkersOnUMAP.pdf",
    width = 18,
    height = 10)
markers_on_umap
dev.off()
png(file = "Output/Plots/20230609_MarkersOnUMAP.png",
    width = 18,
    height = 12,
    units = 'in',
    res = 600,
    pointsize = 14)
markers_on_umap
dev.off()
