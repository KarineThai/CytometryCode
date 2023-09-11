########################################
#' K menas
#' Karine Thai - 2023/08/18
########################################

########################################
#' 0 - Load required packages
########################################

library(vegan)

########################################
#' 1 - Load expression table
########################################

load('Output/RData/PropMFIMetaData.RData')
sub=PropMFIMetaData[PropMFIMetaData$Cohort=='Untreated'&PropMFIMetaData$Group=='MS',]
MS_KM_cal=cascadeKM(sub[,c(1:170)], inf.gr=2,sup.gr=35,iter=100,criterion="ssi") 
plot(MS_KM_cal, sortg=TRUE)

a=data.frame(MS_KM_cal$partition[,3])
head(PropMFIMetaData)
head(a)
class(a)
colnames(a)='4_groups'
Kmeans_4groups=merge(sub[,c(171:210)],a, by=0,all=T)

head(Kmeans_4groups)
Kmeans_4groups=column_to_rownames(Kmeans_4groups, var='Row.names')
table(Kmeans_4groups[,c(3,41)])
