#Obtain MFI boxplots of all markers of interest per cluster for a given comparison 
GetMFIBoxplots = function(df, markers, ID, clusters, metadata, group) {
  #'
  #' @param df Dataframe with cells as rows and parameters as columns. Needs to contain a column with sample IDs and clusters
  #' @param markers Columns in df corresponding to markers to use to generate boxplots
  #' @param ID Column in df containing sample IDs
  #' @param clusters Column in df containing Cluster annotation per cell
  #' @param metadata Dataframe containing metadata to use for comparisons
  #' @param group Group to use for comparison
  #' 
  #' 
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(dplyr, data.table, tibble, ggplot2, ggpubr,ggbeeswarm)
  
  source("~/Documents/PhD/Coding/R_Functions/Median_quartiles.R")
  
  plots = list()
  #If want to scale the columns that will be used and create a dataframe containing scaled parameters, sample IDs, and flowSOM clusters
  df = cbind(data.frame(scale(df[,markers]),data.frame(ID=df[,ID]),data.frame(Clusters=df[,clusters])))
  #Or not scale
  #df = df[,c(markers,ID,clusters)]
  #df = rename(df, ID = ID, Clusters = clusters)
  #Calculate MFI per cluster per sample
  MFIScale = aggregate(df[,-c(ncol(df)-1,ncol(df))], by=list(ID=df[,"ID"], Clusters=df[,"Clusters"]), FUN=median)

  #Split by group
  List_dfPerCluster = split(as.data.table(MFIScale),by="Clusters")
  for(clus in 1:length(List_dfPerCluster)) {
  }
  #Make rownames Sample IDs for MetaData and Each dataframe in the list
  #If Sample ID is the first column in MEtaData
  
  for (clus in 1:length(List_dfPerCluster)) {
    List_dfPerCluster[[clus]]=data.frame(List_dfPerCluster[[clus]])
  }
  
  #metadata=column_to_rownames(metadata, var = colnames(metadata[1]))
  List_dfPerCluster = lapply(List_dfPerCluster, function(x) column_to_rownames(x, var = colnames(x[1])))
  #Add Meta data to each list
  for (clus in 1:length(List_dfPerCluster)) {
    List_dfPerCluster[[clus]]=merge(List_dfPerCluster[[clus]], metadata, by=0)
    }
  #Change Rownames
  List_dfPerCluster = lapply(List_dfPerCluster, function (x) rename(x,Sample_ID = Row.names))
  #Melt table
  for (clus in 1:length(List_dfPerCluster)) {
    List_dfPerCluster[[clus]]=melt(as.data.table(List_dfPerCluster[[clus]]),measure.vars = markers, variable.name = "Marker", value.name = "MFI")
    List_dfPerCluster[[clus]]=data.frame(List_dfPerCluster[[clus]])
  }
  #Create list of plots
  colors_discrete = c("#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#4DBBD5FF", "#B09C85FF", "#3C5488FF")
  for (clus in 1:length(List_dfPerCluster)) {
    title=names(List_dfPerCluster[clus])
    plots[[clus]] = ggplot(subset(List_dfPerCluster[[clus]], !is.na(List_dfPerCluster[[clus]][,group])), aes_string(x="Marker", y="MFI",fill=group)) +  geom_violin(width=0.9,position = position_dodge(width = 0.8)) + stat_summary(fun=median.quartile, position = position_dodge(width=0.8), geom = 'line') + stat_summary(fun="median",size=1, geom="point",position = position_dodge(width=0.8))+geom_quasirandom(size = 0.5,dodge.width = 0.8, varwidth = TRUE, alpha=0.3)+ stat_compare_means(aes_string(group=group), method="t.test",label="p.signif", hide.ns=T) + theme_bw() + ggtitle(title) +  scale_y_continuous(expand = c(.05, 0, .2, 0)) + scale_fill_manual(values=colors_discrete) + theme(axis.text.x = element_text(angle = 45, hjust=1),axis.title.x=element_blank())+theme(legend.position="none")
  }
  return(plots)

}



