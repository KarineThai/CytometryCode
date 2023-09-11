########################################
#' FlowSOM Clustering of Lymphoid and Myeloid cells 
#' Clustering using only cell type markers and using markers of interest
#' Karine Thai - 2023/07/24
########################################

######################
#' 0 - load required packages
######################

library(FlowSOM)
library(flowCore)
library(Biobase)
library(umap)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tidyverse)

######################
#' 1 - Load dataframe of downsample for lympho or myeloid populations
######################

load("Output/RData/Ds_Lympho.RData")
load("Output/RData/Ds_Myeloid.RData")

#####################
#' 2 - Convert dataframe to flowframe
######################

meta <- data.frame(name=colnames(Ds_Myeloid[,1:23]),
                   desc=paste(colnames(Ds_Myeloid[,1:23]))
)
meta$range <- apply(apply(Ds_Myeloid[,1:23],2,range),2,diff)
meta$minRange <- apply(Ds_Myeloid[,1:23],2,min)
meta$maxRange <- apply(Ds_Myeloid[,1:23],2,max)

ff <- new("flowFrame",exprs=as.matrix(Ds_Myeloid[,1:23]),parameters=AnnotatedDataFrame(meta))

######################
#' 3 - Iterate different flowSOM cluster numbers and compute a list of flowSOM objects
######################

# markers <- c("CD161", "HLADR", "CD3", "CD11c", "CD56", "iTCR",  "CD1c", "CD16", "CD20", "CD14")
markers <- c("CD14", "CD11c", "HLADR", "CD123", "CD141", "CD1c", "CD16", "CD3", "CD20")

#allmarkers <- c('SSC.A', 'CD14', 'CD11c', 'CD161', 'HLADR', 'CD3', 'CD56', 'CD18', 'iTCR', 'CD1c', 'CD16', 'CD20', 'DAP12', 'CD69', 'IL7R', 'B2M')
# allmarkers <- c('SSC.A', 'CD14', 'CD11c', 'S100A8A9', 'HLADR', 'CD3', 'CD56', 'CD123', 'CD18', 'CD141', 'CD1c', 'CD16', 'CD20', 'LYZ', 'DAP12', 'CD69', 'IL7R', 'B2M')

# Give importance to iTCR for lympho
# a=c(1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1)
# class(a)
# rbind(a,allmarkers)

#If want to iterate different cluster numbers
source("Functions/IterationNumberClustersFlowSOM.R")
# List_flowSOM_Lympho = IterationNClusFlowSOM(ff = ff, markers = markers, plotdir = 'Output/FlowSOMmary/20230724_FlowSOMmaryLympho_') #This will create a list of FlowSOMs for different numbers of clusters and generate FlowSOMmary for each
List_flowSOM_Myeloid = IterationNClusFlowSOM(ff = ff, markers = markers, plotdir = 'Output/FlowSOMmary/FlowSOMmaryMyeloid_')

#If only want to try one cluster number

# fsom <- FlowSOM(input = ff,
               scale = T,
               compensate=F,
               transform = F,
               colsToUse = markers,
               seed = 1,
               nClus = 20,
               xdim = 10, ydim = 10)


# save(List_flowSOM_Lympho, file = 'Output/RData/List_flowSOM_Lympho.RData')
# save(List_flowSOM_Myeloid, file = 'Output/RData/List_flowSOM_Myeloid.RData')

save(fsom, file = 'Output/RData/flowSOM_Lympho_allmarkers.RData')
save(fsom, file = 'Output/RData/flowSOM_Myeloid_allmarkers.RData')

######################
#' 4 - Add cluster numbers to downsample dataframe and create table containing percent cells in clusters
######################

# LymphoClusRaw=GetMetaclusters(List_flowSOM_Lympho[[4]], meta=NULL)
# Ds_Lympho$LymphoClusRaw=LymphoClusRaw
# LymphoClusPercentages=GetPercentages(List_flowSOM_Lympho[[4]], level='metaclusters')
# write.csv(x = data.frame(LymphoClusPercentages*100), file = 'Output/Tables/PercentLymphoFlowSOM20Clus.csv')

MyeloidClusRaw=GetMetaclusters(List_flowSOM_Myeloid[[4]], meta=NULL)
Ds_Myeloid$MyeloidClusRaw=MyeloidClusRaw
MyeloidClusPercentages=GetPercentages(List_flowSOM_Myeloid[[4]], level='metaclusters')
write.csv(x = data.frame(MyeloidClusPercentages*100), file = 'Output/Tables/PercentMyeloidFlowSOM20Clus.csv')

# ClusterAllMarkers=GetMetaclusters(fsom, meta=NULL)
# Ds_Myeloid$ClusterAllMarkers=ClusterAllMarkers
# ClusPercentages=GetPercentages(fsom, level='metaclusters')
# 
# write.csv(x = data.frame(ClusPercentages*100), file = 'Output/Tables/PercentLymphoFlowSOMAllMarkers.csv')
# write.csv(x = data.frame(ClusPercentages*100), file = 'Output/Tables/PercentMyeloidFlowSOMAllMarkers.csv')

######################
#' 5 - Visualize markers on UMAP and heatmaps
######################

ds = sample_n(Ds_Myeloid, 2000)

#Scale data
scaled_ds = scale(ds[,1:23])
Scaled_Ds_Lympho=scale(Ds_Lympho[,1:23])
Scaled_Ds_Myeloid=scale(Ds_Myeloid[,1:23])

#Markers to use
allmarkersumap <- allmarkers
markersnoB2M <- allmarkers[1:17]

#Generate UMAP

set.seed(1)
UMAP_Lympho = umap(scaled_ds[,allmarkers],n_components=2,  n_neighbors=100,  min_dist=0.4)
UMAP_Myeloid= umap(Scaled_Ds_Myeloid[,allmarkers],n_components=2,  n_neighbors=25,  min_dist=0.4)

# save(UMAP_Lympho, file = 'Output/RData/UMAP_Lympho.RData')
save(UMAP_Myeloid, file = 'Output/RData/UMAP_Myeloid.RData')

# save(UMAP_Lympho, file = 'Output/RData/UMAP_Lympho_allmarkers.RData')
save(UMAP_Myeloid, file = 'Output/RData/UMAP_Myeloid_allmarkers.RData')

#Add UMAP coordinates to Ds dataframe

ds = data.frame(cbind(ds, UMAP_Myeloid$layout))
colnames(ds)[c(31,32)]=c("UMAP1_Allmarkers","UMAP2_Allmarkers")

Ds_Lympho = data.frame(cbind(Ds_Lympho, UMAP_Lympho$layout))
colnames(Ds_Lympho)[c(31,32)]=c("UMAP1_AllMarkers","UMAP2_AllMarkers")
Ds_Myeloid = data.frame(cbind(Ds_Myeloid, UMAP_Myeloid$layout))
colnames(Ds_Myeloid)[c(31,32)]=c("UMAP1_AllMarkers","UMAP2_AllMarkers")

save(Ds_Lympho, file = "Output/RData/Ds_Lympho.RData")
save(Ds_Myeloid, file = "Output/RData/Ds_Myeloid.RData")

#Visualize markers on UMAP

plots_markersclus=list()
for(mark in colnames(ds)[c(1,3,5:23)]){
  plots_markersclus[[mark]]=ggplot(ds,aes_string('UMAP1_Allmarkers','UMAP2_Allmarkers',col=mark))+geom_point(size=0.25)+theme_minimal(base_size = 8.5)+theme(text=element_text(size=14))+ scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(9,"RdYlBu")))
}


plots_markersclus=list()
for(mark in colnames(Ds_Lympho)[c(1,3,5:23)]){
plots_markersclus[[mark]]=ggplot(Ds_Lympho,aes_string('UMAP1_AllMarkers','UMAP2_AllMarkers',col=mark))+geom_point(size=0.25)+theme_minimal(base_size = 8.5)+theme(text=element_text(size=14))+ scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(9,"RdYlBu")))
}
for(mark in colnames(Ds_Myeloid)[c(1,3,5:23)]){
  plots_markersclus[[mark]]=ggplot(Ds_Myeloid,aes_string('UMAP1_AllMarkers','UMAP2_AllMarkers',col=mark))+geom_point(size=0.25)+theme_minimal(base_size = 8.5)+theme(text=element_text(size=14))+ scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(9,"RdYlBu")))
}

plots_allmarkers=list()
allmarkers=colnames(Ds_Myeloid)[c(5:23)]
for(mark in allmarkers){
  plots_allmarkers[[mark]]=ggplot(Ds_Myeloid,aes_string('UMAP1','UMAP2',col=mark))+geom_point(size=0.25)+theme_minimal(base_size = 8.5)+theme(text=element_text(size=14))+ scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(9,"RdYlBu")))
}
#Save Plot
png(file = "Output/Plots/MarkersOnUMAPAllMarkers_Myeloid.png",
    width = 18,
    height = 12,
    units = 'in',
    res = 600,
    pointsize = 14)
cowplot::plot_grid(plotlist = plots_markersclus, nrow = 5)
dev.off()

#Visualize flowSOM clusters on UMAP
load('Output/RData/ColorPaletteP40.RData')
aggreg = aggregate(Ds_Lympho[,c("UMAP1", "UMAP2")], by=list(MyeloidClusRaw=Ds_Myeloid$MyeloidClusRaw2),median)
FlowSOM_Plot=ggplot(Ds_Myeloid,aes(UMAP1,UMAP2))+geom_point(aes(col=MyeloidClusRaw2),size=0.2)+theme_bw()+geom_label_repel(data=aggreg, aes(label=MyeloidClusRaw2, color=MyeloidClusRaw2, fontface='bold'), show.legend = F,max.overlaps = Inf)+scale_color_manual(values = P40)+theme(legend.position = 'none')+ylab('UMAP 2')+xlab('UMAP 1')+theme(axis.title =element_text(size = 16)) 
png(file = "Output/Plots/FlowSOM_Myeloid_UMAP2.png",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 14,
    units = 'in',
    res = 600,
    pointsize = 16)
FlowSOM_Plot
dev.off()

#Visualize markers expression in each cluster by heatmap
MFIs <- GetMetaclusterMFIs(fsom)
Heatmap_Clusters=pheatmap(MFIs[,allmarkers],scale = "column", color = inferno(20), fontsize = 20)
png(file = "Output/Plots/FlowSOMAllMarkers_Myeloid_Heatmap.png",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 8,
    units = 'in',
    res = 600,
    pointsize = 24)
Heatmap_Clusters
dev.off()

#Visualize individual clusters on UMAP
source('Functions/PopFlowSOM_UMAPgrey.R')
PopFlowSOMGrey = ClusFlowSOM_UMAPBW(x = Ds_Lympho, n = 25, PopCol = 'ClusterAllMarkers', UMAP1 = 'UMAP1_AllMarkers', UMAP2='UMAP2_AllMarkers')
png(file = "Output/Plots/FlowSOMAllMarkers_PopGrey_Lympho.png",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 16,
    units = 'in',
    res = 600,
    pointsize = 24)
cowplot::plot_grid(plotlist = PopFlowSOMGrey)
dev.off()

PopFlowSOMGrey = ClusFlowSOM_UMAPBW(x = Ds_Myeloid, n = 20, PopCol = 'MyeloidClusRaw2', UMAP1 = 'UMAP1', UMAP2='UMAP2')
png(file = "Output/Plots/FlowSOMAllMarkers_PopGrey_Myeloid2.png",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 16,
    units = 'in',
    res = 600,
    pointsize = 24)
cowplot::plot_grid(plotlist = PopFlowSOMGrey)
dev.off()

######################
#' 6 - Generate table of frequencies of clusters in samples and visualize by heatmap
#' This is to check whether some clusters are sample-specific
######################

#Contribution of samples to clusters
#Export excel
TableCells=(table(Ds_Lympho[,c(24,30)])/rowSums(table(Ds_Lympho[,c(24,30)])))*100
write.csv(TableCells, 'Output/Tables/TableCellsSampleInClusterLymphoAllMarkers.csv')
HeatmapPerSample=pheatmap(TableCells, scale = 'column',color = inferno(20))
png(file = "Output/Plots/pHeatmap_ClustersvsSamplesLymphoAllMarkers.png",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 14,
    units = 'in',
    res = 600,
    pointsize = 24)
HeatmapPerSample
dev.off()

TableCells=(table(Ds_Myeloid[,c(24,30)])/rowSums(table(Ds_Myeloid[,c(24,30)])))*100
write.csv(TableCells, 'Output/Tables/TableCellsSampleInClusterMyeloidAllMarkers.csv')
pheatmap(TableCells, scale = 'column')
png(file = "Output/Plots/pHeatmap_ClustersvsSamplesMyeloidAllMarkers.png",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 14,
    units = 'in',
    res = 600,
    pointsize = 24)
pheatmap(TableCells, scale = 'column',color = inferno(20))
dev.off()

######################
#' 7 - Clean and rename clusters
######################

load('Output/RData/List_FlowSOM_Lympho.RData')
List_flowSOM_Lympho[[4]] = UpdateMetaclusters(List_flowSOM_Lympho[[4]], newLabels = c('1' = 'Tcells_1', 
                                                                        '2' = 'Tcells_2', 
                                                                        '3' = 'Exclude_4', 
                                                                        '4' = 'Tcells_1', 
                                                                        '5' = 'Exclude_5',
                                                                        '6' = 'CD56dimNK_1', 
                                                                        '7' = 'Tcells_3', 
                                                                        '8' = 'Tcells_4', 
                                                                        '9' = 'CD56dimNK_1', 
                                                                        '10' = 'iNKT',
                                                                        '11' = 'CD56dimNK_2', 
                                                                        '12' = 'Tcells_4', 
                                                                        '13' = 'Tcells_5', 
                                                                        '14' = 'iNKT', 
                                                                        '15' = 'CD56brightNK',
                                                                        '16' = 'Exclude_6', 
                                                                        '17' = 'Bcells_1', 
                                                                        '18' = 'Bcells_2', 
                                                                        '19' = 'Bcells_1', 
                                                                        '20' = 'Bcells_3'))
save(List_flowSOM_Lympho, file = 'Output/RData/List_flowSOM_Lympho.RData')

load('Output/RData/List_flowSOM_Myeloid.RData')
List_flowSOM_Myeloid[[4]] = UpdateMetaclusters(List_flowSOM_Myeloid[[4]], newLabels = c('1' = 'ConventionalDC_CD1c', 
                                                                        '2' = 'ConventionalDC_CD1c', 
                                                                        '3' = 'ClassicalMono_1', 
                                                                        '4' = 'Exclude_1', 
                                                                        '5' = 'ClassicalMono_2',
                                                                        '6' = 'IntermediateMono_1', 
                                                                        '7' = 'IntermediateMono_2', 
                                                                        '8' = 'ClassicalMono_3', 
                                                                        '9' = 'ClassicalMono_3', 
                                                                        '10' = 'NonClassicalMono_1',
                                                                        '11' = 'IntermediateMono_3', 
                                                                        '12' = 'ClassicalMono_4', 
                                                                        '13' = 'ConventionalDC_CD141', 
                                                                        '14' = 'NonClassicalMono_2', 
                                                                        '15' = 'ConventionalDC_DbNeg',
                                                                        '16' = 'NonClassicalMono_2', 
                                                                        '17' = 'Exclude_2', 
                                                                        '18' = 'PlasmacytoidDC', 
                                                                        '19' = 'Basophils', 
                                                                        '20' = 'Exclude_3'))
save(List_flowSOM_Myeloid, file = 'Output/RData/List_flowSOM_Myeloid.RData')
load('Output/RData/FlowSOM_Lympho_allmarkers.RData')
fsom = UpdateMetaclusters(fsom, newLabels = c('1' = 'Tcells_1', 
                                              '2' = 'Tcells_2', 
                                              '3' = 'Bcells_1', 
                                              '4' = 'Bcells_2', 
                                              '5' = 'Bcells_3',
                                              '6' = 'Tcells_3', 
                                              '7' = 'Bcells_4', 
                                              '8' = 'Tcells_4', 
                                              '9' = 'Tcells_5', 
                                              '10' = 'Bcells_1',
                                              '11' = 'Tcells_5', 
                                              '12' = 'CD56dimNK_4', 
                                              '13' = 'CD56dimNK_1', 
                                              '14' = 'CD56dimNK_2', 
                                              '15' = 'Tcells_2',
                                              '16' = 'Tcells_4', 
                                              '17' = 'Exclude_4', 
                                              '18' = 'CD56dimNK_3', 
                                              '19' = 'CD56dimNK_4', 
                                              '20' = 'iNKT',
                                              '21' = 'Tcells_6', 
                                              '22' = 'Tcells_7', 
                                              '23' = 'iNKT', 
                                              '24' = 'Tcells_5', 
                                              '25' = 'CD56brightNK'))
save(fsom, file = 'Output/RData/FlowSOM_Lympho_allmarkers.RData')

load('Output/RData/FlowSOM_Myeloid_allmarkers.RData')
fsom = UpdateMetaclusters(fsom, newLabels = c('1' = 'IntermediateMono_1', 
                                              '2' = 'IntermediateMono_2', 
                                              '3' = 'ClassicalMono_1', 
                                              '4' = 'Exclude_1', 
                                              '5' = 'ClassicalMono_2',
                                              '6' = 'ClassicalMono_3', 
                                              '7' = 'NonClassicalMono_1', 
                                              '8' = 'NonClassicalMono_2', 
                                              '9' = 'ClassicalMono_4', 
                                              '10' = 'ClassicalMono_5',
                                              '11' = 'ClassicalMono_6', 
                                              '12' = 'ClassicalMono_5', 
                                              '13' = 'ClassicalMono_7', 
                                              '14' = 'NonClassicalMono_3', 
                                              '15' = 'ClassicalMono_4',
                                              '16' = 'ClassicalMono_8', 
                                              '17' = 'ConventionalDC_CD1c', 
                                              '18' = 'ClassicalMono_5', 
                                              '19' = 'ConventionalDC_CD1c', 
                                              '20' = 'ConventionalDC_CD141',
                                              '21' = 'ConventionalDC_CD1c_2', 
                                              '22' = 'PlasmacytoidDC', 
                                              '23' = 'Basophils', 
                                              '24' = 'Exclude_2', 
                                              '25' = 'Exclude_3'))
save(fsom, file = 'Output/RData/FlowSOM_Myeloid_allmarkers.RData')
######################
#' 7 - Apply flowSOM to all cells and add clusters to expression tables
######################

load("Output/RData/TableExprsLympho.RData")
load("Output/RData/TableExprsMyeloid.RData")

load('Output/RData/List_FlowSOM_Lympho.RData')
load('Output/RData/List_flowSOM_Myeloid.RData')
load('Output/RData/FlowSOM_Lympho_allmarkers.RData')
load('Output/RData/FlowSOM_Myeloid_allmarkers.RData')

meta <- data.frame(name=colnames(TableExprsMyeloid[,1:23]),
                   desc=paste(colnames(TableExprsMyeloid[,1:23]))
)
meta$range <- apply(apply(TableExprsMyeloid[,1:23],2,range),2,diff)
meta$minRange <- apply(TableExprsMyeloid[,1:23],2,min)
meta$maxRange <- apply(TableExprsMyeloid[,1:23],2,max)

head(meta)

ff = new("flowFrame",exprs=as.matrix(TableExprsMyeloid[,1:23]),parameters=AnnotatedDataFrame(meta))

fsom_Myeloid = NewData(List_flowSOM_Myeloid[[4]], ff)
save(fsom_Myeloid, file = 'Output/RData/fSOM_MainMarkers_Myeloid_TotalCells.RData')
fsom_Myeloid_AllMarkers = NewData(fsom, ff)
save(fsom_Myeloid_AllMarkers, file = 'Output/RData/fSOM_AllMarkers_Myeloid_TotalCells.RData')
fSOMClus=GetMetaclusters(fsom_Myeloid, meta = NULL)
fSOMClusAllMarkers=GetMetaclusters(fsom_Myeloid_AllMarkers, meta = NULL)
TableExprsMyeloid$fSOM_MainMarkers = fSOMClus
TableExprsMyeloid$fSOM_AllMarkers = fSOMClusAllMarkers
save(TableExprsMyeloid, file = "Output/RData/TableExprsMyeloid.RData")
