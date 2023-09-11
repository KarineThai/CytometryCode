########################################
#' Color umap by patient features
#' Karine Thai - 2023/09/06
########################################

######################
#' 0 - Load required packages and functions
######################

library(ggplot2)
library(ggdist)
library(umap)


######################
#' 1 - Load tables
######################

load('Output/RData/DsMainMarkers.RData')
MetaData = read.csv('Input_Tables/ClinicalMetaData_Clean.csv', na.strings = 'NA')


DsMainMarkers_MetaData <- DsMainMarkers %>%
  left_join(MetaData, by = c("Sample_ID" = "Code"))

# The resulting merged_df will have sample information for each observation

DsMainMarkers$TimeToNextRelapse <- MetaData$TimeToNextRelapse[match(DsMainMarkers$Sample_ID, MetaData$Code)]
Ds_Myeloid$TimeToNextRelapse <- MetaData$TimeToNextRelapse[match(Ds_Myeloid$Sample_ID, MetaData$Code)]
ggplot(na.omit(slice_sample(Ds_Myeloid,n=537000)),aes(UMAP1,UMAP2))+geom_point(aes(col=TimeToNextRelapse),size=4)+
  theme_bw()+ylab('UMAP 2')+xlab('UMAP 1')+theme(axis.title =element_text(size = 16))+scale_color_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 6),colors=rev(brewer.pal(11,"RdYlBu"))) 
ds <- DsMainMarkers %>%
  group_by(Sample_ID) %>%
  sample_n(size = 200, replace = TRUE) %>%
  ungroup()

ggplot(na.omit(slice_sample(Ds_Myeloid,n=537000)), aes(x = UMAP1, y = UMAP2)) +
  geom_point() +  # Scatter plot of UMAP coordinates
  geom_density_2d(aes(color = TimeToNextRelapse)) +  # Overlay 2D density plot
  scale_fill_gradient(name = "TimeToNextRelapse", low = "blue", high = "red") +  # Color scale
  theme_minimal()








