library(flowCore)
library(reshape2)
table_param_estimateLogicle = function(ff, channels) {
  #' Returns a table containing all parameters of an estimateLogicle transformList
  #' 
  #' @ff flowFrame
  #' @channels Channels to transform
  #'
  #'
  lgcl= estimateLogicle(ff, channels)
  list_mat=list()
  
  for (col in colnames(ff)[c(5:24)]) {
    l=as.list(environment(lgcl@transforms[[col]]@f))
    list_mat[[col]]=matrix(l[2:6])
  }
  #Création matrice
  mat=list_mat[[1]]
  for (par in 2:length(list_mat)){
    mat=cbind(mat, list_mat[[par]])
    rownames(mat)=c("Param", "w", "t", "m", "a")
    mat=as.data.frame(mat)
  }
  df=t(mat)
  rownames(df)=df[,1]
  df=df[,2:5]
  rownames(df)=gsub("_logicleTransform","", rownames(df))
  return(df)
  }


