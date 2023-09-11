########################################
#' Normalize transformed FCS files
#' Generates files that are transformed and normalized for batch effects, exported as FCS and as RData flowset
#' Karine Thai - 2023/06/07, Updated on 2023/07/24
########################################

######################
#' 0 - load required packages
######################
library(CytoNorm)
library(flowCore)
library(tidyverse)
library(flowCore)
library(ggplot2)
library(reshape2)
library(ggridges)
library(dplyr)
library(stringr)

######################
#' 1 - Retrieve FCS files from reference samples and experimental samples
######################

files <- list.files('FCS/Transformed', pattern = "fcs$")

######################
#' 2 - Organize data for CytoNorm
#' Identify batches and subset into training and validation data
#' Identify channels to normalize and to use for clustering
######################

data <- data.frame(File = files,
                   Path = file.path('FCS/Transformed', files),
                   Type = stringr::str_match(files, "(ref|(([a-z]{2})_([0-9]{3,4})))")[,2],
                   Batch = stringr::str_match(files, "([batch]{5})(10|11|[1-9])")[,1],
                   stringAsFactors = FALSE, truncate_max_range=F)

data$Type2=c(rep(0,192))
for(T in 1:length(data$Type)){
  if (data$Type[T] == "ref"){
    data$Type2[T]="Train"}
  else{
    data$Type2[T] = "Validation"}
}

train_data <- dplyr::filter(data, Type2 == "Train")
validation_data <- dplyr::filter(data, Type2 == "Validation")

#Channels to normalize
ff <- read.FCS(data$Path[1], truncate_max_range=F)
channels <- flowCore::colnames(ff)[c(5:24)]

#Channels to use for FlowSOM clustering
colsForFlowSOM = flowCore::colnames(ff)[c(5,6,10,11,13,15,16,17,18,19)]

######################
#' 3 - Test whether clustering is appropriate
######################

fsom <- prepareFlowSOM(train_data$Path,
                       channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = 10,
                                             scale = FALSE),
                       transformList = NULL,
                       seed = 1)

cvs <- testCV(fsom,
              cluster_values = c(3,5,10,15)) 

cvs$pctgs$`5`
######################
#' 4 - Train the model
######################

model <- CytoNorm.train(files = train_data$Path,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = NULL,
                        FlowSOM.params = list(nCells = 6000, 
                                              xdim = 10,
                                              ydim = 10,
                                              nClus = 5,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE,
                        plot=TRUE)

######################
#' 5 - Normalize the data
######################

CytoNorm.normalize(model = model,
                   files = validation_data$Path,
                   labels = validation_data$Batch,
                   transformList = NULL,
                   transformList.reverse = NULL,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "FCS/Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE,
                   plot=TRUE)

######################
#' 5 - Visualize normalized data by stacked histograms
######################

files_Norm=list.files(path = 'FCS/Normalized',pattern = '.fcs$')
fs_Norm = read.flowSet(files_Norm)

source("Functions/Create_list_expression_table_from_flowset.R")
regex1="(batch)(10|11|[0-9])"
regex2="([a-z]{2,3})_([0-9]{3,4})"
str_match(fs_Norm@phenoData@data$name, "(batch)(10|11|[1-9])")
Tables_exp = Create_list_expression_table_from_flowset(fs_Norm, regex1 = regex1, regex2 = regex2)
names(Tables_exp)

#Verify that all colnames are the same
a=colnames(Tables_exp[[1]])
for (ta in 2:length(Tables_exp)) {
  a=cbind(a,colnames(Tables_exp[[ta]]))
}
all(a==a[,1])

source("Functions/Expression_table_per_marker.R")
Listpermark=Expression_table_per_marker(Tables_exp)

source("Functions/Graph_per_marker_from_list_per_marker.R")
list_plot=Graph_per_marker_from_list_per_marker(Listpermark,nbSD = 8)
names(list_plot)

pdf(file = "Output/Plots/20230608_StackedHistoPerMarkerNorm1.pdf",
    width = 16,
    height = 12) 
cowplot::plot_grid(plotlist=list_plot[c(5:8)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230608_StackedHistoPerMarkerNorm2.pdf",
    width = 16,
    height = 12)
cowplot::plot_grid(plotlist=list_plot[c(9:12)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230608_StackedHistoPerMarkerNorm3.pdf",
    width = 16,
    height = 12)
cowplot::plot_grid(plotlist=list_plot[c(13:16)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230608_StackedHistoPerMarkerNorm4.pdf",
    width = 16,
    height = 12)
cowplot::plot_grid(plotlist=list_plot[c(17:20)], nrow=1)
dev.off()

pdf(file = "Output/Plots/20230608_StackedHistoPerMarkerNorm5.pdf",
    width = 16,
    height = 12)
cowplot::plot_grid(plotlist=list_plot[c(21:24)], nrow=1)
dev.off()
