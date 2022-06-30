library(flowCore)
library(reshape2)
library(ggplot2)
library(ggridges)
iteration_lgcl=function(ff, param, rangeW=c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5)){
  #' Create a list of plots per value of w in tranformLogicle function
  #' 
  #' @param ff le flowframe utilisé pour tester les différentes itérations
  #' @param param paramètre à transformer à écrire entre guillemets
  #' @param rangeW un vecteur de valeurs à essayer pour le W
  #' 
  #' 
  list_vector_transformed=list()
  for (ra in rangeW) {
  transformation = logicleTransform(w = ra, t=262144, m=4.5, a=0)
  translist <- transformList(param, transformation) 
  transformed <- transform(ff,translist)
  list_vector_transformed[[as.character(ra)]]=data.frame(transformed@exprs)[[gsub("-",".",param)]]
  }
  #Création du dataframe
  df = list_vector_transformed[[1]]
  
  for (ra in 2:length(list_vector_transformed)) {
  df = cbind(df,list_vector_transformed[[ra]])
  }
  colnames(df)=names(list_vector_transformed)
  m=melt(df)
  print(head(m))
  p=ggplot(m, aes(x = value, y =as.character(Var2), fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01,show.legend = FALSE) +
    scale_fill_viridis_c(option = "C") +
    xlab("")+
    ylab("")+
    labs(title=param)
  return(p)
}

#Essai
iteration_lgcl(ff=ff, param="FJComp-BV650-A")
