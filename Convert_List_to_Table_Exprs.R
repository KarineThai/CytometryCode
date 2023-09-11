ListToTableExprs = function(x) {
  #'
  #'@param x List of dataframes per sample containing expression matrices
  #'
  #'
  
  x[[1]]$Sample_ID = names(x[1])
  complete_df = x[[1]]
  
  for (sam in 2:length(x)) {
    x[[sam]]$Sample_ID = names(x[sam])
    complete_df = rbind(complete_df, x[[sam]])
  }
  complete_df$Sample_ID = factor(complete_df$Sample_ID)
  complete_df$Total_Cell_ID = 1:nrow(complete_df)
  return(complete_df)
}
