########################################
#' Split list of expression matrices into lymphoid and myeloid cells based on flowSOM metaclusters
#' Karine Thai - 2023/07/24
########################################

######################
#' 1 - Load list of tables of expression with flowSOM MetaCluster
######################

load("Output/RData/ListTablesExpCut_wfSOMMetaClus.RData")

######################
#' 2 - Create new lists of expression for lympho and myeloid cells
#' Lymphoid metaclusters are 1, 2, 4, 5
#' Myeloid metaclusters are 3, 6, 7
######################

List_Exprs_Lympho = list()
List_Exprs_Myeloid = list()
MetaClus_Lympho = c("1","2", "4","5")
MetaClus_Myeloid = c("3","6","7")

#Subset Lympho
for (sam in 1:length(ListTablesExpCut_wfSOMMetaClus)) {
  List_Exprs_Lympho[[sam]]=ListTablesExpCut_wfSOMMetaClus[[sam]][ListTablesExpCut_wfSOMMetaClus[[sam]]$fSOM_MetaClus %in% MetaClus_Lympho,]
}
names(List_Exprs_Lympho)=names(ListTablesExpCut_wfSOMMetaClus)
#Subset Myeloid
for (sam in 1:length(ListTablesExpCut_wfSOMMetaClus)) {
  List_Exprs_Myeloid[[sam]]=ListTablesExpCut_wfSOMMetaClus[[sam]][ListTablesExpCut_wfSOMMetaClus[[sam]]$fSOM_MetaClus %in% MetaClus_Myeloid,]
}
names(List_Exprs_Myeloid)=names(ListTablesExpCut_wfSOMMetaClus)
#Reset rownames
for (i in 1:length(List_Exprs_Lympho)) {
  rownames(List_Exprs_Lympho[[i]]) = seq(length=nrow(List_Exprs_Lympho[[i]]))
  print(head(rownames(List_Exprs_Lympho[[i]])))
}
for (i in 1:length(List_Exprs_Myeloid)) {
  rownames(List_Exprs_Myeloid[[i]]) = seq(length=nrow(List_Exprs_Myeloid[[i]]))
  print(head(rownames(List_Exprs_Myeloid[[i]])))
}

#Make sure that number of rows for Lympho+Myeloid = Total PBMCs
for (sam in 1:length(ListTablesExpCut_wfSOMMetaClus)) {
  print(nrow(List_Exprs_Lympho[[sam]])+nrow(List_Exprs_Myeloid[[sam]])==nrow(ListTablesExpCut_wfSOMMetaClus[[sam]]))
}

save(List_Exprs_Lympho, file = "Output/RData/List_Exprs_Lympho.RData")
save(List_Exprs_Myeloid, file = "Output/RData/List_Exprs_Myeloid.RData")

######################
#' 3 - Convert lists to tables
######################

source("Functions/Convert_List_to_Table_Exprs.R")
TableExprsLympho = ListToTableExprs(List_Exprs_Lympho)
TableExprsMyeloid = ListToTableExprs(List_Exprs_Myeloid)

save(TableExprsLympho, file = "Output/RData/TableExprsLympho.RData")
save(TableExprsMyeloid, file = "Output/RData/TableExprsMyeloid.RData")

######################
#' 4 - Create a downsample of each table
#' Will be used to train flowSOM and UMAP visualizations
######################

source("Functions/TableExprsDownsampleNoCellID.R")
load("Output/RData/List_Exprs_Lympho.RData")
load("Output/RData/List_Exprs_Myeloid.RData")

set.seed(1)
Ds_Lympho=TableExprsDownsampleNoCellID(List_Exprs_Lympho, 3000)
set.seed(1)
Ds_Myeloid=TableExprsDownsampleNoCellID(List_Exprs_Myeloid, 3000)

save(Ds_Lympho, file = "Output/RData/Ds_Lympho.RData")
save(Ds_Myeloid, file = "Output/RData/Ds_Myeloid.RData")

