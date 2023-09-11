#Utilise library(flowCore)
#Fonction de Rose-Marie

Create_list_expression_table_from_flowset=function(x,regex1,regex2){
  #' Create a liste of plots per marker from a files with FCS in
  #' 
  #' @param x le flowframe
  #' @param regex1 regex qui permet d'identifer le sample (partie 1)
  #' @param regex2 regex qui permet d'identifer le sample (partie 2)
  #' 
  #' 

  
  #Creer les listes necessaires pour la fonction
  tables_exp=list() #liste des tables d'expression
  names=c() #liste des antigènes pour chacun des canaux
  codesample=c() #code pour le sample 
  
  
  for(fi in 1:length(x)){
    codesample[fi]=paste0(stringr::str_match(x[[fi]]@description$`$FIL`,regex1)[,1],
                          "_",
                          stringr::str_match(x[[fi]]@description$`$FIL`,regex2)[,1])
    ff=x[[fi]]
    for (par in 1:length(as.vector(ff@parameters@data$desc))){
      if (is.na(as.vector(ff@parameters@data$desc)[par])){
        names[[par]]=as.vector(ff@parameters@data$name)[par]
      }
      else{
        names[[par]]=as.vector(ff@parameters@data$desc)[par] 
      }
    }
    namesV=as.vector(unlist(names))
    newmat=data.frame(ff@exprs)
    colnames(newmat)=namesV
    #newmat$`FSC-A`=NULL
    #newmat$`FSC-H`=NULL
    #newmat$`SSC-A`=NULL
    #newmat$`SSC-H`=NULL
    newmat$`FJComp-IRD800-A`=NULL
    newmat$Time=NULL
    tables_exp[[fi]]=newmat
  }
  names(tables_exp)=codesample
  return (tables_exp)
}
 