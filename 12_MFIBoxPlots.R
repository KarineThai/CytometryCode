########################################
#' Boxplots of MFIs
#' Karine Thai - 2023/08/18
########################################

########################################
#' 0 - Load required packages
########################################

library(dplyr)
library(textshape)
library(data.table)

########################################
#' 1 - Load expression table
########################################

load('Output/RData/TableExprsPBMCMain.RData')
MetaData = read.csv('Input_Tables/ClinicalMetaData_Clean.csv', na.strings = 'NA')
MetaData = textshape::column_to_rownames(MetaData, loc = 1)

MetaData$NEDA5yrs = factor(MetaData$NEDA5yrs, levels = c('NEDA', 'EDA'))
MetaData$Severity5yrs = factor(MetaData$Severity5yrs, levels = c('NEDA', 'MEDA', 'Medium', 'High'))
MetaData$RadiologicalActivity5yrs = factor(MetaData$RadiologicalActivity5yrs, levels = c('Inactive', 'Active'))


source("Functions/BoxPlotsMFI_perMarker_perCluster.R")
mark = c("S100A8A9", 'IL7R')
mark2 = c('LYZ', 'CD69')
mark3 = c('CD18', 'B2M', 'DAP12')
MetaDataF=MetaData[MetaData$Sex=='F',]
MetaDataM=MetaData[MetaData$Sex=='M',]

Plot_list_Group=GetMFIBoxplots(df = TableExprsPBMCMain[TableExprsPBMCMain$fSOM_MetaClus=='3'|TableExprsPBMCMain$fSOM_MetaClus=='6'|TableExprsPBMCMain$fSOM_MetaClus=='7',], markers = mark, ID = "Sample_ID", clusters = "fSOM_MainMarkers", metadata = MetaData, group = "Group")
Plot_list_Progression=GetMFIBoxplots(df = TableExprsPBMCMain, markers = mark, ID = "Sample_ID", clusters = "fSOM_MainMarkers", metadata = MetaData, group = "Sever")
Plot_list_NEDA=GetMFIBoxplots(df = TableExprsPBMCMain[TableExprsPBMCMain$fSOM_MetaClus=='3'|TableExprsPBMCMain$fSOM_MetaClus=='6'|TableExprsPBMCMain$fSOM_MetaClus=='7',], markers = mark, ID = "Sample_ID", clusters = "fSOM_MainMarkers", metadata = MetaData, group = "NEDA5yrs")
Plot_list_Severity=GetMFIBoxplots(df = TableExprsPBMCMain, markers = mark, ID = "Sample_ID", clusters = "fSOM_MainMarkers", metadata = MetaData, group = "Severity5yrs")
Plot_list_Radiological=GetMFIBoxplots(df = TableExprsPBMCMain, markers = mark, ID = "Sample_ID", clusters = "fSOM_MainMarkers", metadata = MetaData, group = "RadiologicalActivity5yrs")
Plot_list_NEDAM=Plot_list_NEDA

ggarrange(plotlist = Plot_list_Severity[15:26], common.legend = T)
