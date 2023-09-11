########################################
#' Cluster visualizations and identity with MEM
#' Karine Thai - 2023/08/16
########################################

######################
#' 0 - Load required packages
######################

library(dplyr)
library(MEM)
library(umap)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(textshape)

######################
#' 1 - Load table of MFIs and proportions with clinical data
######################

load('Output/RData/TableExprsLympho.RData')
load('Output/RData/TableExprsMyeloid.RData')

######################
#' 2 - Bind tables
######################

TableExprsPBMC = rbind(TableExprsLympho, TableExprsMyeloid)
TableExprsPBMC = data.frame(TableExprsPBMC)
TableExprsPBMCMain = TableExprsPBMC[!grepl('Exclude', TableExprsPBMC$fSOM_MainMarkers),]
TableExprsPBMCAll = TableExprsPBMC[!grepl('Exclude', TableExprsPBMC$fSOM_AllMarkers),]
print(unique(TableExprsPBMCMain$fSOM_MainMarkers))
print(unique(TableExprsPBMC$fSOM_AllMarkers))

######################
#' 3 - Add column for cluster number
######################

TableExprsPBMCMain$cluster = TableExprsPBMCMain$fSOM_MainMarkers
TableExprsPBMCAll$cluster = TableExprsPBMCAll$fSOM_AllMarkers
TableExprsPBMCMain$cluster <- as.numeric(factor(TableExprsPBMCMain$cluster, levels = c('Bcells_1',
                                                                                       'Bcells_2',
                                                                                       'Bcells_3',
                                                                                      'Tcells_1',
                                                                                       'Tcells_2',
                                                                                       'Tcells_3',
                                                                                       'Tcells_4',
                                                                                       'Tcells_5',
                                                                                       'iNKT',
                                                                                       'CD56brightNK',
                                                                                       'CD56dimNK_1',
                                                                                       'CD56dimNK_2',
                                                                                       'ClassicalMono_1',
                                                                                       'ClassicalMono_2',
                                                                                       'ClassicalMono_3',
                                                                                       'ClassicalMono_4',
                                                                                       'IntermediateMono_1',
                                                                                       'IntermediateMono_2',
                                                                                       'IntermediateMono_3',
                                                                                       'NonClassicalMono_1',
                                                                                       'NonClassicalMono_2',
                                                                                       'PlasmacytoidDC',
                                                                                       'ConventionalDC_CD1c',
                                                                                       'ConventionalDC_CD141',
                                                                                      'ConventionalDC_DbNeg',
                                                                                       'Basophils')))
TableExprsPBMCAll$cluster <- as.numeric(factor(TableExprsPBMCAll$cluster, levels = c('Bcells_1',
                                                                                     'Bcells_2',
                                                                                     'Bcells_3',
                                                                                     'Bcells_4',
                                                                                     'Tcells_1',
                                                                                     'Tcells_2',
                                                                                     'Tcells_3',
                                                                                     'Tcells_4',
                                                                                     'Tcells_5',
                                                                                     'Tcells_6',
                                                                                     'Tcells_7',
                                                                                     'iNKT',
                                                                                     'CD56brightNK',
                                                                                     'CD56dimNK_1',
                                                                                     'CD56dimNK_2',
                                                                                     'CD56dimNK_3',
                                                                                     'CD56dimNK_4',
                                                                                     'ClassicalMono_1',
                                                                                     'ClassicalMono_2',
                                                                                     'ClassicalMono_3',
                                                                                     'ClassicalMono_4',
                                                                                     'ClassicalMono_5',
                                                                                     'ClassicalMono_6',
                                                                                     'ClassicalMono_7',
                                                                                     'ClassicalMono_8',
                                                                                     'IntermediateMono_1',
                                                                                     'IntermediateMono_2',  
                                                                                     'NonClassicalMono_1',
                                                                                     'NonClassicalMono_2',
                                                                                     'NonClassicalMono_3',
                                                                                     'PlasmacytoidDC',
                                                                                     'ConventionalDC_CD1c',
                                                                                     'ConventionalDC_CD1c_2',
                                                                                     'ConventionalDC_CD141',
                                                                                     'Basophils')))

######################
#' 4 - Downsample for umap
######################
set.seed(123)
DsMainMarkers = sample_n(TableExprsPBMCMain, size = 100000)
set.seed(123)
DsAllMarkers = sample_n(TableExprsPBMCAll, size = 100000)

######################
#' 5 - Compute MEM
######################

mem.res.main = MEM(exp_data = DsMainMarkers,transform = FALSE, markers='5,6,8,9,10,11,12,14,15,16,17,18', zero.ref = FALSE)
build.heatmaps(
  mem.res.main, 
  cluster.MEM = "both",
  cluster.medians = "none",
  display.thresh = 1,
  newWindow.heatmaps = FALSE,
  output.files = TRUE,
  labels = TRUE,
  only.MEMheatmap = TRUE,
  output.dir = 'Output/MEM'
)

mem.res.all = MEM(exp_data = DsAllMarkers,transform = FALSE, markers='5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23', zero.ref = FALSE)
build.heatmaps(
  mem.res.all, 
  cluster.MEM = "both",
  cluster.medians = "none",
  display.thresh = 1,
  newWindow.heatmaps = FALSE,
  output.files = TRUE,
  labels = TRUE,
  only.MEMheatmap = TRUE,
  output.dir = 'Output/MEM'
)
TableExprs$ClusterLabel <- ClusterLabelsDF$ClusterLabel[match(TableExprs$cluster_number_column, ClusterLabelsDF$ClusterNumber)]

######################
#' 6 - Add MEM to table Exprs
######################
load('Output/RData/TableExprsPBMCMain.RData')

ClusterMEMsDF = read.csv('Input_Tables/fSOM_MainMarkers_MatchedMEM.csv')
ClusterMEMsDF$MEM = gsub(' ', '', ClusterMEMsDF$MEM)
TableExprsPBMCMain$MEM <- ClusterMEMsDF$MEM[match(TableExprsPBMCMain$fSOM_MainMarkers, ClusterMEMsDF$fSOM_MainMarkers)]
TableExprsPBMCMain$Cluster_MEM = paste(TableExprsPBMCMain$fSOM_MainMarkers, TableExprsPBMCMain$MEM, sep = ' | ')

save(TableExprsPBMCMain, file = 'Output/TableExprs_Pheno_Clean.RData')


######################
#' 7 - Scale data for UMAP
######################

Scaled_DsMainMarkers=scale(DsMainMarkers[,1:23])

######################
#' 8 - Generate UMAP with scaled parameters
######################

Scaled_UMAPMainMarkers = umap(Scaled_DsMainMarkers[,c(1,3,5,6,8:12,14:18)],n_components=2,  n_neighbors=15,  min_dist=0.6)

######################
#' 9 - Add UMAP coordinates to table
######################

DsMainMarkers = data.frame(cbind(DsMainMarkers, Scaled_UMAPMainMarkers$layout))
colnames(DsMainMarkers)[c(30,31)]=c("UMAP1","UMAP2")
save(DsMainMarkers, file = 'Output/RData/DsMainMarkers.RData')
save(DsAllMarkers, file = 'Output/RData/DsAllMarkers.RData')

######################
#' 10 - Plot and save clusters on UMAP
######################
load('Output/RData/ColorPaletteP40.RData')
aggreg = aggregate(DsMainMarkers[,c("UMAP1", "UMAP2")], by=list(cluster=DsMainMarkers$cluster),median)
FlowSOM_Plot=ggplot(DsMainMarkers,aes(UMAP1,UMAP2))+geom_point(aes(col=as.factor(cluster)),size=0.5)+theme_bw()+geom_label_repel(data=aggreg, aes(label=cluster, color=as.factor(cluster), fontface='bold'),size=7, show.legend = F,max.overlaps = Inf)+scale_color_manual(values = sample(P40, 26))+theme(legend.position = 'none')+ylab('UMAP 2')+xlab('UMAP 1')+theme(axis.title =element_text(size = 16)) 
png(file = "Output/Plots/FlowSOM_PBMC_UMAP_clusternumber.png",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 14,
    units = 'in',
    res = 600,
    pointsize = 16)
FlowSOM_Plot
dev.off()

######################
#' 11 - Plot and save markers on UMAP
######################

molecules=colnames(DsMainMarkers)[c(3,5:23)]

plots=list()
for(mol in molecules){
  plots[[mol]]=ggplot(DsMainMarkers,aes_string("UMAP1","UMAP2",col=mol))+geom_point(size=0.25)+theme_minimal(base_size = 8.5)+theme(text=element_text(size=14))+ scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(9,"RdYlBu")))
}

markers_on_umap = cowplot::plot_grid(plotlist=plots[c(1,7,11,5,15,8,14,2,3,6,13,12,9,4,10,16:20)], nrow=4)

png(file = "Output/Plots/MarkersOnUMAP_PBMCMainMarkers.png",
    width = 18,
    height = 12,
    units = 'in',
    res = 600,
    pointsize = 14)
markers_on_umap
dev.off()

######################
#' 12 - Visualize individual clusters on UMAP
######################

source('Functions/PopFlowSOM_UMAPgrey.R')
PopFlowSOMGrey = ClusFlowSOM_UMAPBW(x = DsMainMarkers, n = 26, PopCol = 'cluster', UMAP1 = 'UMAP1', UMAP2='UMAP2')
png(file = "Output/Plots/FlowSOMMainMarkers_PopGrey_PBMC.png",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 16,
    units = 'in',
    res = 600,
    pointsize = 24)
cowplot::plot_grid(plotlist = PopFlowSOMGrey)
dev.off()

######################
#' 13 - Heatmap of markers per cluster
######################

MFI = aggregate(TableExprsPBMCMain[,c(5:23)], by=list(Clusters=TableExprsPBMCMain[,"fSOM_MainMarkers"]), FUN=median)
MFI = textshape::column_to_rownames(MFI, loc = 1)
Heatmap_Clusters = pheatmap(MFI,color = blues9, fontsize = 20)
png(file = "Output/Plots/FlowSOMMainMarkers_PBMC_Heatmap_Clusternames.png",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10,
    units = 'in',
    res = 600,
    pointsize = 24)
Heatmap_Clusters
dev.off()


