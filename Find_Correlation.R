find_cor=function(df,feature, col){
  #'
  #'@param df df to use with proportions
  #'@param feature column containing clinical factor to correlate with proportions
  #'@param col vector of columns corresponding to clusters
  #'
  #'
  
  col_cluster = colnames(df)[col]
  est = c()
  pvals = c()
  for(col in col_cluster){
    res = cor.test(df[,col],df[,feature])
    est=c(est,res$estimate)
    pvals=c(pvals,res$p.value)
  }
  out = data.frame(cbind(est,pvals,col_cluster))
  out[,1]=as.numeric(out[,1])
  out[,2]=as.numeric(out[,2])
  print(head(out[order(out[,'pvals']),]))
  out[,'col_cluster']=factor(out[,'col_cluster'],out[order(-out[,'pvals']),"col_cluster"])
  return(out)
}
