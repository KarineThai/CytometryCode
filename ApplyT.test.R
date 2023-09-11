ApplyT.test=function(df,feature, col){
  #'
  #'@param df df to use with proportions
  #'@param feature column containing two clinical groups
  #'@param col vector of columns corresponding to clusters
  #'
  #'
  
  variable = colnames(df)[col]
  log2fc=c()
  pvals = c()
  for(col in variable){
    res = t.test(as.numeric(df[,col])~df[,feature], data = df, var.equal=T)
    log2fc=c(log2fc,log2(res[["estimate"]][[names(res[["estimate"]])[2]]]/res[["estimate"]][[names(res[["estimate"]])[1]]]))
    pvals=c(pvals,res$p.value)
  }
  out = data.frame(cbind(log2fc,pvals,variable))
  out[,1]=as.numeric(out[,1])
  out[,2]=as.numeric(out[,2])
  out[,'p.adj']=p.adjust(out[,'pvals'], method='fdr', n=length(out[,'pvals']))
  out[,4]=as.numeric(out[,4])
  print(head(out[order(out[,'p.adj']),]))
  out[,'variable']=factor(out[,'variable'],out[order(-out[,'p.adj']),"variable"])
  return(out)
}
