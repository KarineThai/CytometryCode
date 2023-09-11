library(CytoNorm)
library(flowCore)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggridges)
library(ggcyto)
library(xlsx)

#Create flowSet for transformed files not normalized and normalized+transformed files
FCSLocation_Norm = '/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/Normalized_FCS'
fs_Norm = flowCore::read.flowSet(path = FCSLocation_Norm, pattern="fcs$", truncate_max_range=T)

#Subset cells for which values for each parameters are 0 to 10
GateList = list()
channels=colnames(fs_Norm[[1]])[5:24]

for (i in 1:length(channels )){
  para<-channels[i]
  GateList[[i]]<- rectangleGate(par=c(0,5))
  GateList[[i]]@parameters[[1]]@parameters=para
}
res=filter(fs_Norm, GateList[[1]]&GateList[[2]]&GateList[[3]]&GateList[[4]]&GateList[[5]]&GateList[[6]]&GateList[[7]]&GateList[[8]]&GateList[[9]]&GateList[[10]]&GateList[[11]]&GateList[[12]]&GateList[[13]]&GateList[[14]]&GateList[[15]]&GateList[[16]]&GateList[[17]]&GateList[[18]]&GateList[[19]]&GateList[[20]])

fs_Norm_Filtered=Subset(fs_Norm,res)
#Rename FCS with wrong label
fs_Norm_Filtered@frames$Norm_batch4_mspatient_stain_gb_mg_128_20140812.fcs@description$`$FIL`="batch4_mspatient_stain_gb_mg_202_20140812.fcs"
fs_Norm_Filtered@frames$Norm_batch4_mspatient_stain_jl_mg_202_20140508.fcs@description$`$FIL`="Norm_batch4_mspatient_stain_jl_mg_185_20140508.fcs"
fs_Norm_Filtered@frames$Norm_batch4_hlctrl_stain_cd_159_20210503.fcs@description$`$FIL`="Norm_batch4_hlctrl_stain_cp_159_20210503.fcs"

save(fs_Norm_Filtered, file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/RData/fs_Norm_Filtered.RData")

setwd("/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/")

source("/Users/p0088821/Documents/PhD/Coding/R_Functions/Create_list_expression_table_from_flowset.R")

#Identify the regex
regex1="(batch)(10|11|[0-9])"
regex2="([a-z]{2,3})_([0-9]{3,4})"
stringr::str_match(fs_Norm_Filtered@phenoData@data$name, "(batch)(10|11|[1-9])")
Tables_exp=Create_list_expression_table_from_flowset(fs_Norm_Filtered,regex1=regex1, regex2=regex2 )
names(Tables_exp)

#Verify that all colenames is the same
a=colnames(Tables_exp[[1]])

for (ta in 2:length(Tables_exp)) {
  a=cbind(a,colnames(Tables_exp[[ta]]))
}

#Verify that all colunm in a is identical
all(a==a[,1])


source("Functions/Expression_table_per_marker.R")

Listpermark=Expression_table_per_marker(Tables_exp)


source("Functions/Graph_per_marker_from_list_per_marker.R")
list_plot=Graph_per_marker_from_list_per_marker(Listpermark, nbSD=8)
names(list_plot)
pdf(file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/Plots/20230608_StackedHistoPerMarkerNormFilter1.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches

cowplot::plot_grid(plotlist=list_plot[c(5:8)], nrow=1)

dev.off()

pdf(file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/Plots/20230608_StackedHistoPerMarkerNormFilter2.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches

cowplot::plot_grid(plotlist=list_plot[c(9:12)], nrow=1)

dev.off()

pdf(file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/Plots/20230608_StackedHistoPerMarkerNormFilter3.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches

cowplot::plot_grid(plotlist=list_plot[c(13:16)], nrow=1)

dev.off()

pdf(file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/Plots/20230608_StackedHistoPerMarkerNormFilter4.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches

cowplot::plot_grid(plotlist=list_plot[c(17:20)], nrow=1)

dev.off()

pdf(file = "/Users/p0088821/Documents/PhD/MSBiobank_Biomarker/Integration_181samples/Output/Plots/20230608_StackedHistoPerMarkerNormFilter5.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches

cowplot::plot_grid(plotlist=list_plot[c(21:24)], nrow=1)

dev.off()
