TableExprsDownsampleNoCellID = function(x, nevents) {
  #'
  #'@param x List of dataframes per sample containing expression matrices
  #'@param nevents The number of events to be randomly subsetted for each sample
  #'
  if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)
  }
  
  x[[1]]$Sample_ID = names(x[1])
  complete_df = sample_n(x[[1]], nevents)
  
  for (sam in 2:length(x)) {
    x[[sam]]$Sample_ID = names(x[sam])
    complete_df = rbind(complete_df, sample_n(x[[sam]], nevents))
  }
  complete_df$Sample_ID = factor(complete_df$Sample_ID)
  return(complete_df)
}
