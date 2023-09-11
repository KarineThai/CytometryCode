#Faire une fonction pour savoir le nombre de cellules dans une liste de table d'expression

NumberCellsPerSampleTruncated=function(x_before,x_after){

  #' @param  x_before, la liste de table donnés pas la fonction Create_list_expression_table_from_flowset
  #' @param  x_after, la liste 

  Vec_names=c()
  Vec_total_cell_before=c()
  Vec_total_cell_after=c()


  Vec_names=names(x_after)
  
  for (nam in Vec_names){
    Vec_total_cell_before[nam]=nrow(data.frame(x_before[nam]))
    Vec_total_cell_after[nam]=nrow(data.frame(x_after[nam]))
  }
  df=data.frame("tot_cell_before"=Vec_total_cell_before,"tot_cell_after"=Vec_total_cell_after, "diff_cell"=(Vec_total_cell_before-Vec_total_cell_after),"Percent_cut"=(100*(Vec_total_cell_before-Vec_total_cell_after)/Vec_total_cell_before))
  options(scipen = 999)
  Total=c(sum(Vec_total_cell_before), 
          sum(Vec_total_cell_after), 
          sum((Vec_total_cell_before-Vec_total_cell_after)),
          100*(sum(Vec_total_cell_before)-sum(Vec_total_cell_after))/sum(Vec_total_cell_before))
  df=rbind(df,Total=Total)
  
  return(df)
}
