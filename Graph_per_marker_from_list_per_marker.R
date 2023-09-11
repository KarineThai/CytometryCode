library(flowCore)
library(ggplot2) 
library(reshape2)
library(ggridges)


Graph_per_marker_from_list_per_marker=function(x, nbSD){
  #' Create a liste of plots per marker from a list given by the fonction "Expression_table_per_marker"
  #' 
  #' @param x list de table d'expression tel que donné par la fonction "Expression_table_per_marker"
  #' @param nbSD nombre d'ecart type à montrer sur le graph 
  #' 
  
  plot_list=list()
  for (ta in 1:length(x)){
  tt_melt=melt(x[[ta]])
  sdtt=sd(tt_melt$value,na.rm=TRUE)
  moytt=mean(tt_melt$value,na.rm=TRUE)
  tt_melt_subset=subset(tt_melt, value>=(moytt-nbSD*sdtt) & value<(moytt+nbSD*sdtt))
  p=ggplot(tt_melt_subset, aes(x = value, y = variable, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01,show.legend = FALSE) +
    scale_fill_viridis_c(option = "C") + 
    xlab("")+
    ylab("")+
    labs(title=names(x)[[ta]])
  #je suis rendu la ajouter mon plot dans ma liste de plot
  plot_list[[names(x)[[ta]]]]=p
  }
return (plot_list)
}
