########################################
#' Marker Enrichment Modeling for FlowSOM clusters
#' Karine Thai - 2023/07/27
########################################

######################
#' 0 - Load required packages
######################

library(MEM)
library(dplyr)

######################
#' 1 - Load dataframe with cells as rows and the last column as cluster number
######################

load('Output/RData/Ds_Lympho.RData')
load('Output/RData/Ds_Myeloid.RData')

######################
#' 2 - Format dataframe
#' Cluster Number needs to be the last column of the dataframe, in numeric format
######################

Ds_Lympho$LymphoClusRaw=as.numeric(as.vector(Ds_Lympho$LymphoClusRaw))
Ds_Lympho <- Ds_Lympho %>% relocate(LymphoClusRaw, .after = UMAP2)
Ds_Lympho$cluster = Ds_Lympho$LymphoClusRaw

Ds_Myeloid$MyeloidClusRaw=as.numeric(as.vector(Ds_Myeloid$MyeloidClusRaw))
Ds_Myeloid <- Ds_Myeloid %>% relocate(MyeloidClusRaw, .after = UMAP2)
Ds_Myeloid$cluster = Ds_Myeloid$MyeloidClusRaw
######################
#' 3 - Compute marker enrichment modeling
######################

mem.res.lympho = MEM(exp_data = Ds_Lympho,transform = FALSE, markers='5,6,8,9,10,11,14,16,17,18', zero.ref = TRUE)
build.heatmaps(
  mem.res.lympho, 
  cluster.MEM = "both",
  cluster.medians = "none",
  display.thresh = 5,
  newWindow.heatmaps = FALSE,
  output.files = TRUE,
  labels = TRUE,
  only.MEMheatmap = TRUE,
  output.dir = 'Output'
)

mem.res.myeloid = MEM(exp_data = Ds_Myeloid,transform = FALSE, markers='5,6,9,10,12,15,16,17,18', zero.ref = TRUE)
build.heatmaps(
  mem.res.myeloid, 
  cluster.MEM = "both",
  cluster.medians = "none",
  display.thresh = 5,
  newWindow.heatmaps = FALSE,
  output.files = TRUE,
  labels = TRUE,
  only.MEMheatmap = TRUE,
  output.dir = 'Output'
)


######################
#' 3 - Add to table expression and clean table
######################

ClusterMEMsDF = read.csv('Input/Tables/Pheno_ClusterNumberWithMEM.csv')

TableExprs$MEM <- ClusterMEMsDF$MEM[match(TableExprs$cluster, ClusterMEMsDF$ClusterNumber)]

save(TableExprs, file = 'Output/TableExprs_Pheno_Clean.RData')

