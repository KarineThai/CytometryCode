########################################
#' Plot histogram distribution of markers in clusters to check if appropriate to use MFIs
#' Karine Thai - 2023/09/05
########################################

######################
#' 0 - Load required packages and functions
######################

library(data.table)
source("Functions/PlotDensityPerCluster.R")
source("Functions/GetShapiroAllClus.R")
######################
#' 1 - Load table of MFIs and proportions with clinical data
######################

load('Output/RData/TableExprsPBMCMain.RData')

######################
#' 2 - Make histograms
######################

dfMFI = List_dfPerCluster[[1]]
for (clus in 2:length(List_dfPerCluster)) {
  dfMFI=rbind(dfMFI, List_dfPerCluster[[clus]])
}
mark=c("S100A8A9","CD18","LYZ","DAP12","CD69","IL7R","B2M")
PlotsMFIinSamples=PlotDensityPerCluster(df=TableExprsPBMCMain,colCluster="fSOM_MainMarkers",markers=mark, group='Group')
ggarrange(plotlist = PlotsMFIinSamples[1:3], common.legend = T)
PlotsMFIinSamplesProg=PlotDensityPerCluster(df=na.omit(dfMFI),colCluster="Clusters",markers=mark, group="Progression_5yrs")
ggarrange(plotlist = PlotsMFIinSamplesProg, common.legend = T)
PlotsMFIinSamplesRad=PlotDensityPerCluster(df=na.omit(dfMFI),colCluster="Clusters",markers=mark, group="Radiological_activity")
ggarrange(plotlist = PlotsMFIinSamplesRad, common.legend = T)
PlotsMFIinSamplesClin=PlotDensityPerCluster(df=na.omit(dfMFI),colCluster="Clusters",markers=mark, group="Clinical_Activity")
ggarrange(plotlist = PlotsMFIinSamplesClin, common.legend = T)

MFIScale = aggregate(TableExprsPBMCMain[,c(1:23)], by=list(ID=TableExprsPBMCMain[,"Sample_ID"], Clusters=TableExprsPBMCMain[,"fSOM_MainMarkers"]), FUN=median)
TableExprsPBMCMain$Group=TableExprsPBMCMain$Sample_ID
TableExprsPBMCMain$Group <- sub("^(mspatient|hlctrl).*", "\\1", TableExprsPBMCMain$Group)
