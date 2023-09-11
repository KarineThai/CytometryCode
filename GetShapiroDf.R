#Function Shapiro test
GetShapiroDf = function(df, clusters, markers, group, group1,group2) {
  #' @param Listdf Dataframe for one cluster containing MFI per sample
  #' @param clusters Column in dataframe containing cluster
  #' @param markers List of markers to test with Shapiro test
  #' @param group Column containing groups to test
  #' @param group1 Name of first group to test
  #' @param group2 Name of seconde group to test
  #' 
  #' 
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ggplot2, ggpubr, ggridges,tidyverse,dplyr)
  
  List_Shap=list()
  List_pvalues=list()
  for (ma in mark) {
      shap=tapply(df[,ma], df[,group], shapiro.test)
      List_Shap[[ma]]=shap
      List_pvalues[[ma]]=data.frame("Marker"=names(List_Shap[ma]), "p.value.group1"=List_Shap[[ma]][["HC"]][["p.value"]],"p.value.group2"=List_Shap[[ma]][["MS"]][["p.value"]])
      #print(colnames(List_pvalues[[ma]]))
      }

  df_Shap=List_pvalues[[1]]
  for (ma in 2:length(List_pvalues)) {
    df_Shap = rbind(df_Shap,List_pvalues[[ma]])
  }
  colnames(df_Shap) = c("Marker", paste0("p.value.",group1), paste0("p.value.", group2))
  df_Shap[,"Cluster"]=df[1:nrow(df_Shap),clusters]
  #print(df_Shap)
  return(df_Shap)
}
#List_Shap1=GetShapiroDf(df=List_dfPerCluster[[1]],clusters = "Clusters", markers=mark, group="Group", group1="HC", group2="MS")
