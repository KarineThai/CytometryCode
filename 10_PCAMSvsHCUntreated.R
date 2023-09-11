########################################
#' PCA to compare MS vs HC
#' Karine Thai - 2023/08/16
########################################

######################
#' 0 - Load required packages
######################

library(tidyr)
library(dplyr)
library(ggfortify)
library(ggforce)
library(factoextra)
library(corrplot)
source("Functions/Median_quartiles.R")

######################
#' 1 - Load table of MFIs and proportions with clinical data
######################

load('Output/RData/PropMFIMetaData.RData')
load('Output/RData/PropTableAll.RData')

######################
#' 2 - Replace NAs in MS_Type column for HC and Untreated for Untreated
#' Clean to remove clusters to exclude
######################

PropTableAll['MSTypeSampling'][is.na(PropTableAll['MSTypeSampling'])] = 'HC'
PropMFIMetaData['MSTypeSampling'][is.na(PropMFIMetaData['MSTypeSampling'])] = 'HC'

PropTableAll <- PropTableAll %>%
  mutate(TreatmentSampling = ifelse(Group == "MS" & is.na(TreatmentSampling), "Untreated", TreatmentSampling))
PropMFIMetaData <- PropMFIMetaData %>%
  mutate(TreatmentSampling = ifelse(Group == "MS" & is.na(TreatmentSampling), "Untreated", TreatmentSampling))

PropTableAll <- PropTableAll %>% select(-contains("Exclude"))
PropMFIMetaData <- PropMFIMetaData %>% select(-contains("Exclude"))

PropTableAll$LastRelapseToSampling=as.Date(PropTableAll$LastRelapseToSampling)
PropTableAll$Date_Sampling=as.Date(PropTableAll$Date_Sampling)
PropTableAll$TimeSinceLastRelapse=as.numeric(difftime(PropTableAll$Date_Sampling, PropTableAll$LastRelapseToSampling, units = 'days'))/30
PropTableAll <- PropTableAll %>% relocate(TimeSinceLastRelapse, .after = LastRelapseToSampling)
PropTableAll['mspatient_mg_108', 'TimeToSPMS']=5.04383561643836
PropTableAll['mspatient_pd_400', 'TimeToSPMS']=2.47945205479452
PropMFIMetaData$LastRelapseToSampling=as.Date(PropMFIMetaData$LastRelapseToSampling)
PropMFIMetaData$Date_Sampling=as.Date(PropMFIMetaData$Date_Sampling)
PropMFIMetaData$TimeSinceLastRelapse=as.numeric(difftime(PropMFIMetaData$Date_Sampling, PropMFIMetaData$LastRelapseToSampling, units = 'days'))/30
PropMFIMetaData <- PropMFIMetaData %>% relocate(TimeSinceLastRelapse, .after = LastRelapseToSampling)
PropMFIMetaData['mspatient_mg_108', 'TimeToSPMS']=5.04383561643836
PropMFIMetaData['mspatient_pd_400', 'TimeToSPMS']=2.47945205479452
save(PropTableAll, file = 'Output/RData/PropTableAll.RData')
save(PropMFIMetaData, file = 'Output/RData/PropMFIMetaData.RData')

######################
#' 3 - Subset untreated cohort with all 
######################

PropTableAllUntreated = PropTableAll[PropTableAll$Cohort == 'Untreated'|PropTableAll$Group=='HC',]
PropMFIMetaDataUntreated = PropMFIMetaData[PropMFIMetaData$Cohort == 'Untreated'|PropTableAll$Group=='HC',]

PropMFIMSOutcome = PropMFIMetaData[PropMFIMetaData$Cohort=='Outcome5yrs'&PropMFIMetaData$Group=='MS',]
######################
#' 4 - Compute PCA
######################
Untreated = PropMFIMetaData[PropMFIMetaData$Group=='MS'&PropMFIMetaData$TreatmentSampling=='Untreated'&PropMFIMetaData$Cohort=='Outcome5yrs',]
PCAMSvsHC = prcomp(PropTableAllUntreated[,c(1:35)], center = T, scale. = T)
PCAPropAll = prcomp (PropMFIMetaDataUntreated[,c(1:170)], center = T, scale. = T)
summary(PCAMSvsHC)
summary (PCAPropAll)
autoplot(PCAPropAll, x=1,y=2, data = PropMFIMetaDataUntreated, color = 'Group') #+ geom_text(aes(label = rownames(PropTableAllUntreated)))

var <- get_pca_var(PCAPropAll, )
corrplot(var$cos2[,c(1:10)], is.corr=FALSE)
PCAType = prcomp(PropMFIMetaData[PropMFIMetaData$Cohort=='Untreated',c(1:170)], scale. = T, center =T)
autoplot(PCAType, x=1,y=3,data = PropMFIMetaData[PropMFIMetaData$Cohort=='Untreated',], color = 'MSTypeSampling')+geom_point(cex=2, aes(color=MSTypeSampling))+ scale_color_brewer(palette = 'Dark2')+theme_bw()+theme(text = element_text(size = 24), axis.text = element_text(size=20)) 

PCAMS = prcomp(PropMFIMSOutcome[,c(1:170)], scale. = T, center = T)
PCAGA = prcomp(PropMFIMetaData[PropMFIMetaData$Group=='MS'&PropMFIMetaData$TreatmentSampling=='Glatiramer Acetate'&PropMFIMetaData$Cohort=='Outcome5yrs',c(1:170)], center = T, scale. = T)
PCAUntreated = prcomp(Untreated[,c(1:170)], center = T, scale. = T)
PCAIFNb = prcomp(PropMFIMetaData[PropMFIMetaData$Group=='MS'&PropMFIMetaData$TreatmentSampling=='IFNb'&PropMFIMetaData$Cohort=='Outcome5yrs',c(1:170)], center = T, scale. = T)
PCAFingo = prcomp(PropMFIMetaData[PropMFIMetaData$Group=='MS'&PropMFIMetaData$TreatmentSampling=='Fingolimod'&PropMFIMetaData$Cohort=='Outcome5yrs',c(1:170)], center = T, scale. = T)
png(filename = 'Output/Plots/PCAMSOutcome_ColorTreatment.png',
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
autoplot(PCAMS, data = PropMFIMSOutcome, color = 'TreatmentSampling')+geom_point(cex=2, aes(color=TreatmentSampling))+ scale_color_npg()+theme_bw()+theme(text = element_text(size = 24), axis.text = element_text(size=20)) 
dev.off()

######################
#' 5 - Make EDSS categories
######################

PropMFIMetaData$EDSS_Sampling_Category=PropMFIMetaData$EDSS_Sampling
PropMFIMetaData$EDSS_Sampling_Category=gsub("1$|1.5$|2$","1-2",PropMFIMetaData$EDSS_Sampling_Category)
PropMFIMetaData$EDSS_Sampling_Category=gsub("2.5$|3$|3.5$|4$|4.5$|6$","2+",PropMFIMetaData$EDSS_Sampling_Category)
Untreated$PC1= PCAUntreated$x[rownames(Untreated),1]
Untreated$PC2 = PCAUntreated$x[rownames(Untreated),2]
colors_discrete = c("#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#4DBBD5FF", "#B09C85FF", "#3C5488FF")
my_comparisons <- list( c("0", "1-2"), c("0", "2+"), c("1-2", "2+") )
ggplot(Untreated, aes(x=as.character(EDSS_Sampling_Category), y=PC1, fill=EDSS_Sampling_Category)) +  geom_violin(width=0.9,position = position_dodge(width = 0.8)) +stat_summary(fun=median.quartile, position = position_dodge(width=0.8), geom = 'line') + stat_summary(fun="median",size=1, geom="point",position = position_dodge(width=0.8))+ geom_quasirandom(size = 0.5,dodge.width = 0.7, varwidth = TRUE, alpha=0.3) + stat_compare_means(comparisons = my_comparisons, method="t.test",label="p.signif") + theme_bw() + scale_fill_manual(values=colors_discrete)




summary(PCAGA)
autoplot(PCAUntreated, x=2,y=1,data = PropMFIMetaData[PropMFIMetaData$Group=='MS'&PropMFIMetaData$TreatmentSampling=='Untreated'&PropMFIMetaData$Cohort=='Outcome5yrs',], color = 'dEDSS_5yrs')+theme_bw()
