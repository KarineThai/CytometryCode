########################################
#' FlowSOM Clustering of all PBMCs using only general cell population markers
#' Generate 7 clusters to divide lymphoid and mmyeloid cells
#' Karine Thai - 2023/07/24
########################################

######################
#' 0 - load required packages
######################

library(FlowSOM)
library(Biobase)
library(flowCore)
library(ggplot2)
library(viridis)
library(umap)
library(pheatmap)
library(data.table)

######################
#' 1 - Load list of tables of expression
######################

load("Output/RData/ListTablesExpCut.RData")

######################
#' 2 - Create a downsampled table using 3000 cells per sample
######################

source("Functions/TableExprsDownsampleNoCellID.R")
source("Functions/Create_list_expression_table_from_flowset.R")

set.seed(1)
Ds_PBMC_3000cellsPerSample=TableExprsDownsampleNoCellID(ListTablesExpCut, 3000)
# save(Ds_PBMC_3000cellsPerSample, file = 'Output/RData/Ds_PBMC_3000cellsPerSample.RData')

######################
#' 3 - Convert Ds_PBMC_3000cellsPerSample to a flowframe for input to flowSOM
######################

load('Output/RData/Ds_PBMC_3000cellsPerSample.RData')
meta <- data.frame(name=colnames(Ds_PBMC_3000cellsPerSample[,1:23]),
                   desc=paste(colnames(Ds_PBMC_3000cellsPerSample[,1:23]))
)
meta$range <- apply(apply(Ds_PBMC_3000cellsPerSample[,1:23],2,range),2,diff)
meta$minRange <- apply(Ds_PBMC_3000cellsPerSample[,1:23],2,min)
meta$maxRange <- apply(Ds_PBMC_3000cellsPerSample[,1:23],2,max)

ff <- new("flowFrame",exprs=as.matrix(Ds_PBMC_3000cellsPerSample[,1:23]),parameters=AnnotatedDataFrame(meta))

######################
#' 4 - Compute the flowSOM object and save flowSOMmary file
#' We use a downsample from all samples to train the flowSOM and will apply the map to all cells
######################

SOM_x <- 5
SOM_y <- 5
markers_of_interest <- c("FSC.A", "SSC.A", "CD14", "CD11c", "HLADR", "CD3", "CD56", "CD123", "CD16", "CD20")

fsom <- FlowSOM(input = ff,
                scale = T,
                compensate=F,
                transform = F,
                colsToUse = markers_of_interest,
                seed = 1,
                nClus = 7,
                xdim = SOM_x, ydim = SOM_y)
FlowSOMmary(fsom, plotFile = "Output/FlowSOMmary/20230725_FlowSOMmary_meta_fsom.pdf")

fSOM_MetaClus_Ds=GetMetaclusters(fsom, meta=NULL)
Ds_PBMC_3000cellsPerSample$fSOM_MetaClus_Ds=fSOM_MetaClus_Ds
save(Ds_PBMC_3000cellsPerSample, file = 'Output/RData/Ds_PBMC_3000cellsPerSample.RData')

######################
#' 5 - Apply flowSOM to all cells
######################

#Convert list of tables of expression to a dataframe with sample IDs and absolute cell IDs
source('Functions/Convert_List_to_Table_Exprs.R')
TableExprsPBMC=ListToTableExprs(ListTablesExpCut)
save(TableExprsPBMC, file = 'Output/RData/TableExprsPBMC.RData')

#Convert dataframe to flowframe
metaPBMC <- data.frame(name=colnames(TableExprsPBMC[,1:23]),
                       desc=paste(colnames(TableExprsPBMC[,1:23]))
)
metaPBMC$range <- apply(apply(TableExprsPBMC[,1:23],2,range),2,diff)
metaPBMC$minRange <- apply(TableExprsPBMC[,1:23],2,min)
metaPBMC$maxRange <- apply(TableExprsPBMC[,1:23],2,max)

ffPBMC <- new("flowFrame",exprs=as.matrix(TableExprsPBMC[,1:23]),parameters=AnnotatedDataFrame(metaPBMC))

#Apply flowSOM and save
meta_fsom = NewData(fsom, ffPBMC)
save(meta_fsom, file = "Output/RData/Meta_fsom_AllCells.RData")
FlowSOMmary(meta_fsom, plotFile = "Output/FlowSOMmary/20230725_FlowSOMmary_meta_fsom.pdf")

######################
#' 6 - Plot outliers
######################

outlier_report <- TestOutliers(meta_fsom,
                               madAllowed = 4,
                               channels = meta_fsom$map$colsUsed,
                               plotFile = "Output/Plots/20230627_Outlier_meta_fsom.pdf")

######################
#' 7 - Add cluster numbers to table expression
######################

fSOM_MetaClus=GetMetaclusters(meta_fsom, meta=NULL)
TableExprsPBMC$fSOM_MetaClus=fSOM_MetaClus
save(TableExprsPBMC, file = 'Output/RData/TableExprsPBMC.RData')

######################
#' 8 - Split back dataframe into list of dataframes by sample, for subsequent steps
######################

ListTablesExpCut_wfSOMMetaClus = split(as.data.table(TableExprsPBMC), by='Sample_ID')
save(ListTablesExpCut_wfSOMMetaClus, file="Output/RData/ListTablesExpCut_wfSOMMetaClus.RData")

######################
#' 9 - Visualizations of clusters and markers by heatmap
######################

load('Output/RData/Meta_fsom_AllCells.RData')
mfis <- GetMetaclusterMFIs(meta_fsom)
colnames(mfis)
Pheatmap_fSOMMetaClus_AllMarkers=pheatmap(mfis[,c(1,3,5:23)],scale = "column", color = inferno(20), fontsize = 14)
Pheatmap_fSOMMetaClus_MetaMarkers=pheatmap(mfis[,c(1,3,5,6,9,10,11,12,17,18)],scale = "column", color = inferno(20), fontsize = 14)

png(file = "Output/Plots/Heatmap_FlowSOM_MetaClus_AllMarkers.png",
    width = 8,
    height = 5,
    units = 'in',
    res = 600,
    pointsize = 24)
Pheatmap_fSOMMetaClus_AllMarkers
dev.off()
png(file = "Output/Plots/Heatmap_FlowSOM_MetaClus_MetaMarkers.png",
    width = 6,
    height = 4,
    units = 'in',
    res = 600,
    pointsize = 24)
Pheatmap_fSOMMetaClus_MetaMarkers
dev.off()

######################
#' 10 - Visualize clusters and markers on UMAP
#' Make new downsample because training flowSOM was not saved and clusters could not be traced back to cells in initial downsample
#' Compute UMAP using only markers used to cluster
######################

load('Output/RData/ListTablesExpCut_wfSOMMetaClus.RData')

#Downsample
set.seed(1)
Ds_PBMC_3000cellsPerSample = TableExprsDownsampleNoCellID(x = ListTablesExpCut_wfSOMMetaClus, nevents = 3000)

#Scale
Scaled_Ds_PBMC_3000cellsPerSample = scale(Ds_PBMC_3000cellsPerSample[,1:23])

#Generate UMAP
set.seed(1)
UMAP_PBMC_MetaClus = umap(Scaled_Ds_PBMC_3000cellsPerSample[,markers_of_interest],n_components=2,  n_neighbors=15,  min_dist=0.6)
save(UMAP_PBMC_MetaClus, file = 'Output/RData/UMAP_PBMC_MetaClus.RData')

#Add UMAP coordinates to Ds_PBMC_3000cellsPerSample
Ds_PBMC_3000cellsPerSample = data.frame(cbind(Ds_PBMC_3000cellsPerSample, UMAP_PBMC_MetaClus$layout))
colnames(Ds_PBMC_3000cellsPerSample)[c(26,27)]=c("UMAP1","UMAP2")
#save(Ds_PBMC_3000cellsPerSample, file = 'Output/RData/Ds_PBMC_3000cellsPerSample.RData')

#Visualize individual clusters on UMAP
source('Functions/PopFlowSOM_UMAPgrey.R')
PopFlowSOMGrey = ClusFlowSOM_UMAPBW(x = Ds_PBMC_3000cellsPerSample, n = 7, PopCol = 'fSOM_MetaClus_Ds', UMAP1 = 'UMAP1', UMAP2='UMAP2')
png(file = "Output/20230723_FlowSOM_MetaClus_PopGrey.png",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10,
    units = 'in',
    res = 600,
    pointsize = 24)
cowplot::plot_grid(plotlist = PopFlowSOMGrey)
dev.off()

#Visualize markers on UMAP
plots=list() 
for (mark in markers_of_interest) {
  plots[[mark]] = ggplot(Ds_PBMC_3000cellsPerSample, aes_string('UMAP1', 'UMAP2', col = mark)) + geom_point(size=0.25) + theme_minimal(base_size = 8.5) + theme(text=element_text(size=14))+ scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(9,"RdYlBu")))
}

#Save Plot
png(file = "Output/Plots/20230723_MarkersOnUMAP_MetaClusMarkers.png",
    width = 14,
    height = 8,
    units = 'in',
    res = 600,
    pointsize = 14)
cowplot::plot_grid(plotlist = plots, nrow = 3)
dev.off()
