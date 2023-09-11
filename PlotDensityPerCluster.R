#Fonction pour voir la distribution de chaque marqueur par cluster en histogrammes afin de déterminer si on peut faire le MFI
PlotDensityPerCluster = function(df, colCluster, markers, group) {
  #'@param df Dataframe containing expression values per marker
  #'@param colCluster Column in df containing the cluster numbers
  #'@param markers markers to plot
  #'@param group Column containing the group that will be used to color the curves
  
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(ggplot2, ggpubr, ggridges,tidyverse,dplyr)
  
  plots = list()
  for (ma in markers) {
    p = ggplot(df, aes_string(x=ma, y=colCluster, color=group)) + aes(fill=stat(x)) +  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +scale_fill_viridis_c(option = "C")  + theme(axis.title.y = element_blank())+ ylab("")+scale_colour_manual(values=c("black","white"))+theme_bw()
    plots[[ma]]=p 
     }
  return(plots)
}
