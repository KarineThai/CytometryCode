Expression_table_per_marker=function(x,nevents=20000){
  #' Create a liste of plots per marker from am  table create by Create_list_expression_table
  #' 
  #' @param x list de table d'expression tel que donné par la fonction "create_list_expresion_table"
  #' @param nevents  nombre d'evenement pour la creation des plots
  #' 
  list_m_e=list()
  
  a=colnames(x[[1]])
  for (mark in 1:length(a)){
    #print(mark)
    tt=x[[1]][a[[mark]]][1:nevents,]
    for (ta in 2:length(x)){
      #print(ta)
      tt=cbind(data.frame(tt),x[[ta]][a[[mark]]][1:nevents,])
    }
    colnames(tt)=names(x)
    list_m_e[[a[[mark]]]]=tt 
  }
  return (list_m_e)
}  


