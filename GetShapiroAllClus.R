GetShapiroAllClus=function(Listdf){
  #' @param Listdf List of df of MFI per sample per Cluster
  #' 
  #' 
  
  source("~/Documents/PhD/Coding/R_Functions/GetShapiroDf.R")
  
  ListShapiroPerCluster = list()
  for (clus in 1:length(List_dfPerCluster)) {
    ListShapiroPerCluster[[clus]]=GetShapiroDf(df=List_dfPerCluster[[clus]],clusters = "Clusters", markers=mark, group="Group", group1="HC", group2="MS")
  }
  df_Shap=ListShapiroPerCluster[[1]]
  for (clus in 2:length(ListShapiroPerCluster)) {
    df_Shap = rbind(df_Shap,ListShapiroPerCluster[[clus]])
  }
  return(df_Shap)
}

