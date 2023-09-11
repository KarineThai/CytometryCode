library(flowCore)
library(reshape2)

logicle_transformList = function(ff, df, colw) {
  #' Create a transformList using logicle tranformation with custom parameters per channel
  #'
  #' @ff flowFrame
  #' @df a data frame containing the channels as rownames and a column containing the w value associated to each channel
  #' @colw the name of the column containing the w values
  
  #Create a list of parameters and of w values
  par=rownames(df)
  wpar=df[[colw]]
  lgcl=list()
  for (pa in par) {
    lgcl[[pa]] = logicleTransform(w = wpar[[which(par==pa)]], t=262144, m=4.5, a=0)
  }
  transformList = transformList(par, lgcl)
  for (pa in par) {
  print(pa)  
  print(as.list(environment(transformList@transforms[[pa]]@f)))
  }
  return(transformList)
}



