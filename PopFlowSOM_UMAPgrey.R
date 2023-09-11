#Faire la distribution des clusters flowSOM en noir sur un UMAP gris
ClusFlowSOM_UMAPBW=function(x,n, PopCol, UMAP1, UMAP2){
  #' Create a liste of plots per marker from a files with FCS in
  #'
  #' @param x data frame avec une colonne qui contient les différent clusters
  #' @param n nombre de différent clusters
  #' @param PopCol nom de la colonne avec les cluster à mettre en string
  #' @param UMAP1 nom de la colonne comprenant parametre UMAP X1 en string
  #' @param UMAP2 nom de la colonne comprenant parametre UMAP X2 en string
  #'
  plots=list()
  for(pop in 1:n){
    title=paste("Cluster",pop)
    plots[[pop]]=ggplot(x,aes_string(UMAP1,UMAP2,col=factor(x[PopCol]==pop,level=c(TRUE,FALSE))))+geom_point(size=0.25)+theme_bw()+scale_color_grey()+theme(legend.position = "none", plot.title = element_text(size = 20),axis.title.x = element_blank(),
                                                                                                                                                            axis.title.y = element_blank(), 
                                                                                                                                                            axis.text.x = element_blank(),
                                                                                                                                                            axis.text.y = element_blank(),
                                                                                                                                                            axis.ticks = element_blank())+ggtitle(title)
  }
  return(plots)
}
