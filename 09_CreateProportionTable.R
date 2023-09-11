########################################
#' Make proportion of clusters per sample
#' Karine Thai - 2023/08/02
########################################

######################
#' 0 - Load required packages
######################

library(data.table)
library(dplyr)
library(textshape)

######################
#' 1 - Load table expression and table containing groups
######################
load("Output/RData/TableExprsLympho.RData")
load("Output/RData/TableExprsMyeloid.RData")
MetaData = read.csv('Input_Tables/ClinicalMetaData_Clean.csv', na.strings = 'NA')

######################
#' 2 - Create proportion tables
######################

#For Lympho main markers
PropTableLymphoMain = table(TableExprsLympho[,c("fSOM_MainMarkers","Sample_ID")])
for(i in 1:ncol(PropTableLymphoMain)){
  PropTableLymphoMain[,i]=(PropTableLymphoMain[,i]/sum(PropTableLymphoMain[,i]))*100
}
PropTableLymphoMain=data.frame(PropTableLymphoMain)
PropTableLymphoMain=data.frame(dcast(data.table(PropTableLymphoMain),Sample_ID~fSOM_MainMarkers, value.var="Freq"))
PropTableLymphoMain=column_to_rownames(PropTableLymphoMain, loc = 1)

#For Myeloid main markers
PropTableMyeloidMain = table(TableExprsMyeloid[,c("fSOM_MainMarkers","Sample_ID")])
for(i in 1:ncol(PropTableMyeloidMain)){
  PropTableMyeloidMain[,i]=(PropTableMyeloidMain[,i]/sum(PropTableMyeloidMain[,i]))*100
}
PropTableMyeloidMain=data.frame(PropTableMyeloidMain)
PropTableMyeloidMain=data.frame(dcast(data.table(PropTableMyeloidMain),Sample_ID~fSOM_MainMarkers, value.var="Freq"))
PropTableMyeloidMain=column_to_rownames(PropTableMyeloidMain, loc = 1)

#For Lympho all markers
PropTableLymphoAll = table(TableExprsLympho[,c("fSOM_AllMarkers","Sample_ID")])
for(i in 1:ncol(PropTableLymphoAll)){
  PropTableLymphoAll[,i]=(PropTableLymphoAll[,i]/sum(PropTableLymphoAll[,i]))*100
}
PropTableLymphoAll=data.frame(PropTableLymphoAll)
PropTableLymphoAll=data.frame(dcast(data.table(PropTableLymphoAll),Sample_ID~fSOM_AllMarkers, value.var="Freq"))
PropTableLymphoAll=column_to_rownames(PropTableLymphoAll, loc = 1)

#For Myeloid all markers
PropTableMyeloidAll = table(TableExprsMyeloid[,c("fSOM_AllMarkers","Sample_ID")])
for(i in 1:ncol(PropTableMyeloidAll)){
  PropTableMyeloidAll[,i]=(PropTableMyeloidAll[,i]/sum(PropTableMyeloidAll[,i]))*100
}
PropTableMyeloidAll=data.frame(PropTableMyeloidAll)
PropTableMyeloidAll=data.frame(dcast(data.table(PropTableMyeloidAll),Sample_ID~fSOM_AllMarkers, value.var="Freq"))
PropTableMyeloidAll=column_to_rownames(PropTableMyeloidAll, loc = 1)

######################
#' 3 - Format Metadata table
######################

MetaData = column_to_rownames(MetaData, loc = 1)

######################
#' 4 - Add group to prop table
######################

PropTableMain = merge(PropTableLymphoMain, PropTableMyeloidMain, by = 0)
PropTableMain = column_to_rownames(PropTableMain, loc = 1)
PropTableMain = merge(PropTableMain, MetaData, by = 0)
PropTableMain = column_to_rownames(PropTableMain, loc = 1)

PropTableAll = merge(PropTableLymphoAll, PropTableMyeloidAll, by = 0)
PropTableAll = column_to_rownames(PropTableAll, loc = 1)
PropTableAll = merge(PropTableAll, MetaData, by = 0)
PropTableAll = column_to_rownames(PropTableAll, loc = 1)


save(PropTableMain, file = 'Output/RData/PropTableMain.RData')
save(PropTableAll, file = 'Output/RData/PropTableAll.RData')

######################
#' 5 - Get MFIs of markers of interest in each cluster
#' Markers of interest in lympho: CD18, DAP12, CD69, IL7R, B2M
#' Markers of interest in myeloid: S100A8A9, CD18, LYZ, DAP12, CD69, IL7R, B2M
######################

load('Output/RData/TableExprsLympho.RData')
load('Output/RData/TableExprsMyeloid.RData')
# For Lympho
dfLympho = TableExprsLympho[,c('CD18', 'DAP12', 'CD69', 'IL7R', 'B2M', 'Sample_ID', 'fSOM_MainMarkers')] #Table of expression subsetted for markers of interest, cluster name, and sample ID columns only
dfLympho = dfLympho[!grepl("Exclude", dfLympho$fSOM_MainMarkers),] #Remove clusters to exclude
unique(dfLympho$fSOM_MainMarkers) #Make sure that we have all clusters

MFILympho = aggregate(dfLympho[,-c(ncol(dfLympho)-1,ncol(dfLympho))], by=list(Sample_ID=dfLympho[,"Sample_ID"], fSOM_MainMarkers=dfLympho[,"fSOM_MainMarkers"]), FUN=median) #Obtenir les MFI de chaque marqueur par sample et par cluster

MFILympho=data.frame(melt(data.table(MFILympho), id.vars=1:2,measure.vars = 3:7, variable.name='Marker',value.name = 'MFI'))
MFILympho$Marker_in_Cluster=paste(MFILympho[,'Marker'],MFILympho[,'fSOM_MainMarkers'], sep='_')
MFILympho=MFILympho[,c(1,4,5)] #Enlever les colonnes marqueur et cluster puisque l'on vient de créer une nouvelle colonne qui merge les deux
MFILympho=reshape(MFILympho, idvar = 'Sample_ID', timevar = 'Marker_in_Cluster', direction='wide') #Obtenir une longue table où chaque colonne est le MFI dechaque marqueur dans chaque cluster
MFILympho=column_to_rownames(MFILympho, loc = 1)

# For Myeloid
dfMyeloid = TableExprsMyeloid[,c('S100A8A9', 'CD18', 'LYZ', 'DAP12', 'CD69', 'IL7R', 'B2M', 'Sample_ID', 'fSOM_MainMarkers')] #Table of expression subsetted for markers of interest, cluster name, and sample ID columns only
dfMyeloid = dfMyeloid[!grepl("Exclude", dfMyeloid$fSOM_MainMarkers),] #Remove clusters to exclude
unique(dfMyeloid$fSOM_MainMarkers) #Make sure that we have all clusters

MFIMyeloid = aggregate(dfMyeloid[,-c(ncol(dfMyeloid)-1,ncol(dfMyeloid))], by=list(Sample_ID=dfMyeloid[,"Sample_ID"], fSOM_MainMarkers=dfMyeloid[,"fSOM_MainMarkers"]), FUN=median) #Obtenir les MFI de chaque marqueur par sample et par cluster

MFIMyeloid=data.frame(melt(data.table(MFIMyeloid), id.vars=1:2,measure.vars = 3:8, variable.name='Marker',value.name = 'MFI'))
MFIMyeloid$Marker_in_Cluster=paste(MFIMyeloid[,'Marker'],MFIMyeloid[,'fSOM_MainMarkers'], sep='_')
MFIMyeloid=MFIMyeloid[,c(1,4,5)] #Enlever les colonnes marqueur et cluster puisque l'on vient de créer une nouvelle colonne qui merge les deux
MFIMyeloid=reshape(MFIMyeloid, idvar = 'Sample_ID', timevar = 'Marker_in_Cluster', direction='wide') #Obtenir une longue table où chaque colonne est le MFI dechaque marqueur dans chaque cluster
MFIMyeloid=column_to_rownames(MFIMyeloid, loc = 1)

######################
#' 6 - Merge MFI tables
######################
MFIinClusters = merge(MFILympho, MFIMyeloid, by = 0)
MFIinClusters = column_to_rownames(MFIinClusters, loc = 1)
head(MFIinClusters)

######################
#' 7 - Merge MFI tables
######################

load('Output/RData/PropTableMain.RData')

PropMFIMetaData = merge(PropTableMain, MFIinClusters, by = 0)
PropMFIMetaData = column_to_rownames(PropMFIMetaData, loc = 1)
PropMFIMetaData = PropMFIMetaData %>% relocate(c(33:72), .after = MFI.IL7R_Basophils)

save(PropMFIMetaData, file = 'Output/RData/PropMFIMetaData.RData')
