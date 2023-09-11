AddMetaClusToListTableExpression = function(x, fSOMObject,regex1,regex2) {
  #' Add Metaclusters and Cell ID to expression table for all FlowFrames
  #' 
  #' @param x list of table of expression (fct Create_list_expression_table_from_flowset)
  #' @param fSOMObject flowSOM Object containing metaclusters for all cells in all samples
  #' @param regex1 regex qui permet d'identifer le sample (partie 1)
  #' @param regex2 regex qui permet d'identifer le sample (partie 2)
  #' The combination of regex1_regex2 needs to be the same as the sample codes in x 
  #' 
#Get MetaCluster list 
  List_MetaCluster = GetMetaclusters(fSOMObject)
  List_Clusters_Per_Sample = list()
  codesample=c() 
  
  for(fi in 1:length(fSOMObject$metaData)){
    codesample[fi]=paste0(stringr::str_match(names(fSOMObject$metaData[fi]),regex1)[,1],
                          "_",
                          stringr::str_match(names(fSOMObject$metaData[fi]),regex2)[,1])
    List_Clusters_Per_Sample[[fi]]= List_MetaCluster[fSOMObject$metaData[[fi]][1]:fSOMObject$metaData[[fi]][2]]
    names(List_Clusters_Per_Sample)=codesample
    x[[names(List_Clusters_Per_Sample[fi])]]$Cell_ID=fSOMObject$metaData[[fi]][1]:fSOMObject$metaData[[fi]][2]
  }
  for (nam in 1:length(List_Clusters_Per_Sample)) {
    print(dim(x[[names(List_Clusters_Per_Sample[nam])]]))
    print(dim(as.data.frame(List_Clusters_Per_Sample[[nam]])))
    print(names(List_Clusters_Per_Sample[nam]))
    x[[names(List_Clusters_Per_Sample[nam])]]$fSOM_MetaClus=List_Clusters_Per_Sample[[nam]]
  }
  return(x)
}
  
  