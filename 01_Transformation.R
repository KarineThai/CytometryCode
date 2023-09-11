########################################
#' Transform FCS files to a linear scale
#' Generate a flowset containing transformed files
#' Karine Thai - 2023/06/07, Updated on 2023/07/24
########################################

######################
#' 0 - Load required packages
######################

library(flowCore)
library(reshape2)
library(readxl)
library(tidyverse)
library(stringr)

######################
#' 1 - Import tables containing transformation parameters
#' Path is the directory where the table containing w parameters for logicle transformation is stored
######################

Path='Input_Tables/Table_Parameters_EstimateLogicle.xlsx'
ParamTable = read_excel(Path, col_names =T)
ParamTable=column_to_rownames(ParamTable,var='...1')

######################
#' 2 - Read FCS files as flowset
#' FCSLocation is the directory where FCS files are stored
######################

files = list.files(FCSLocation, pattern='.fcs$')
fs = read.flowSet(files=files, path = FCSLocation) 

######################
#' 3 - Create transformList with optimal values for logicleTransform
######################

source("Functions/Custom_TransformList.R")
ff = fs[[1]]
transList = logicle_transformList(ff, ParamTable, "w_iteration")

######################
#' 4 - Apply transformList to flowSet
#' Save flowSet as R object and as transformed FCS files
######################

fs_Transformed = transform(fs, transList)
save(fs_Transformed, file = "Output/RData/fs_Transformed.RData")
write.flowSet(fs_Transformed, outdir = "FCS/Transformed")

######################
#' 4 - Visualize transformations by stacked histograms
######################

#Create a list of expression matrices by sample from flowset
source("Functions/Create_list_expression_table_from_flowset.R")
regex1="(batch)(10|11|[0-9])"
regex2="([a-z]{2,3})_([0-9]{3,4})"
str_match(fs@phenoData@data$name, "(batch)(10|11|[1-9])")
Tables_exp=Create_list_expression_table_from_flowset(fs_Transformed,regex1=regex1, regex2=regex2 )

#Verify that all colnames are the same
a=colnames(Tables_exp[[1]])
for (ta in 2:length(Tables_exp)) {
  a=cbind(a,colnames(Tables_exp[[ta]]))
}
all(a==a[,1])

#Create expression table per marker
source("Functions/Expression_table_per_marker.R")
Listpermark=Expression_table_per_marker(Tables_exp)

#Create a graph of stacked histograms of samples for each marker
source("Functions/Graph_per_marker_from_list_per_marker.R")
list_plot=Graph_per_marker_from_list_per_marker(Listpermark,nbSD = 8)

names(list_plot)
pdf(file = "Output/Plots/20230607_StackedHistoPerMarkerTransformed1.pdf", 
    width = 16, 
    height = 12) 
cowplot::plot_grid(plotlist=list_plot[c(5:8)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230607_StackedHistoPerMarkerTransformed2.pdf",   
    width = 16, 
    height = 12) 
cowplot::plot_grid(plotlist=list_plot[c(9:12)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230607_StackedHistoPerMarkerTransformed3.pdf", 
    width = 16,
    height = 12)
cowplot::plot_grid(plotlist=list_plot[c(13:16)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230607_StackedHistoPerMarkerTransformed4.pdf",   
    width = 16,
    height = 12)
cowplot::plot_grid(plotlist=list_plot[c(17:20)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230607_StackedHistoPerMarkerTransformed5.pdf",
    width = 16,
    height = 12)
cowplot::plot_grid(plotlist=list_plot[c(21:24)], nrow=1)
dev.off()


