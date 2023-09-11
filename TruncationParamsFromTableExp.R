TruncationParamsFromTableExp <- function(x, table_min_max,min, max, marker) {
 
  #'Cut outliers form table_expression based on min and max values 
  #' @param x list of table_expression
  #' @param table_min_max dataframe with min and max values for each marker
  #' @param min name of the colunm with min values
  #' @param max name of the colunm with max values
  #' @param marker name of the colunm with marker names
  #' @return list of table_expression without outliers
library(dplyr) 

  rownames(table_min_max) <- table_min_max[,marker]
  for (mkr in rownames(table_min_max)){
    x <- lapply(x, function(ta) ta %>%
                  filter(!!sym(mkr) >= table_min_max[mkr,min] & !!sym(mkr) <= table_min_max[mkr,max]))
    print(head(x[[ta]]))
    } 
return(x)
}
