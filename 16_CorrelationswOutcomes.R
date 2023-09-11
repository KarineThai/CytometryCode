########################################
#' T test of cluster and MFI changes
#' Karine Thai - 2023/08/16
########################################

######################
#' 0 - Load required packages and functions
######################

library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggbeeswarm)
library(RColorBrewer)
source("Functions/ApplyT.test.R")
source('Functions/makeStars.R')
source('Functions/Find_Correlation_Spearman.R')

######################
#' 1 - Load table of MFIs and proportions with clinical data
######################

source("Functions/ApplyT.test.R")
source('Functions/Find_Correlation.R')
load('Output/RData/PropMFIMetaData.RData')
load('Output/RData/PropTableAll.RData')

######################
#' 2 - Update factors and levels
######################

PropMFIMetaData$Progression_5yrs=factor(PropMFIMetaData$Progression_5yrs, levels=c('Stable','Progressor'))
PropMFIMetaData$Group=factor(PropMFIMetaData$Group, levels=c('HC','MS'))
PropMFIMetaData$RadiologicalActivity5yrs=factor(PropMFIMetaData$RadiologicalActivity5yrs, levels=c('Inactive','Active'))
PropMFIMetaData$ClinicalActivity5yrs=factor(PropMFIMetaData$ClinicalActivity5yrs, levels=c('ClinicallyInactive','ClinicallyActive'))
PropMFIMetaData$NEDA5yrs=factor(PropMFIMetaData$NEDA5yrs, levels=c('NEDA','EDA'))
PropMFIMetaData$Severity5yrs=factor(PropMFIMetaData$Severity5yrs, levels=c('NEDA','MEDA','Medium','High'))


######################
#' 3 - T test
######################
res_group = ApplyT.test(PropMFIMetaData[PropMFIMetaData$TreatmentSampling=='Untreated'|PropMFIMetaData$Group=='HC',], 'Group', c(1:170))
res=ApplyT.test(PropMFIMetaData[PropMFIMetaData$Cohort=='Outcome5yrs'&PropMFIMetaData$Group=='MS',], "NEDA5yrs", c(1:170))
cor_Lesions = find_cor(PropMFIMetaData, 'NbNewLesions5yrs', c(1:170))
cor_ARR = find_cor(PropMFIMetaData, 'ARR_5yrs', c(1:170))

cor_TimeToRelapse_Pop = find_cor(PropTableAll, 'TimeToNextRelapse', c(1:35))
cor_Onset = find_cor(PropMFIMetaData[PropMFIMetaData$Group=='MS'&PropMFIMetaData$TreatmentSampling=='Untreated',], 'TimeSinceOnset', c(1:170))

cor_Pop = find_cor(PropTableAll, 'dEDSS_5yrs', c(1:35))
cor_Age = find_cor(PropMFIMetaData[PropMFIMetaData$Group=='MS'&PropMFIMetaData$TreatmentSampling=='Untreated',], 'Age', c(1:170))
cor_TimeSinceRelapse = find_cor(PropMFIMetaData, 'TimeSinceLastRelapse', c(1:170))
cor_TimeToSPMS = find_cor(PropMFIMetaData, 'TimeToSPMS', c(1:170))

ggplot(PropMFIMetaData,aes(x=ClassicalMono_4,y=TimeToNextRelapse))+geom_point(aes(col=TreatmentSampling))+stat_cor(size=7)+geom_smooth(method='lm')+theme_minimal() +theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18)) + scale_color_npg()
head(cor_TimeToRelapse[order(cor_TimeToRelapse[,'pvals']),])
colnames(PropTableAll)
res = find_cor(PropTableAll, 'TimeSinceOnset', c(1:35))
ggplot(cor_TimeToRelapse[cor_TimeToRelapse$pvals<0.05,],aes(est,reorder(col_cluster,est)))+geom_col(width=.75,col='black',aes(fill=-log10(pvals)))+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')), name="-log10(p value)")+xlab("Correlation (r)")+ylab("")+theme_minimal()+theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16), legend.text = element_text(size=14), legend.title=element_text(size=16))       

######################
#' 4 - Find correlation using spearman
######################

spcor_TTNR = find_cor_spearman(PropMFIMetaData, 'TimeToNextRelapse', c(1:170))
spcr_ARR = find_cor_spearman(PropMFIMetaData, 'Age', c(1:170))

######################
#' 4 - Plots of correlation for time to next relapse
######################

#Extract most significant parameters that predict time to next relapes
cor_TimeToRelapse = find_cor(PropMFIMetaData, 'TimeToNextRelapse', c(1:170))
cor_TimeToRelapse$p.signif = makeStars(cor_TimeToRelapse$pvals)
png(filename = 'Output/Plots/TimeToNextRelapseSignificantParametersCorrelation.png',
    width = 13,
    height = 7,
    units = 'in',
    res = 600)
ggplot(cor_TimeToRelapse[cor_TimeToRelapse$pvals<0.05,],aes(est,reorder(col_cluster,est)))+
  geom_col(width=.75,col='black',aes(fill=-log10(pvals)))+
  scale_fill_gradientn(colors=brewer.pal(name = 'YlOrRd', n = 9), name="-log10(p value)")+
  xlab("Correlation (r) with time to next relapse (years)")+ylab("")+
  geom_text(aes(label = p.signif,  hjust = ifelse(est < 0, 1.1, -0.1)), size=10, vjust = 0.8)+
  theme_bw()+
  theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16), legend.text = element_text(size=14), legend.title=element_text(size=16))       
dev.off()

#Make correlation plots of most significant parameters with Time to next relapse (iNKT, IntermediateMono_2, ClassicalMono_1, ClassicalMono_4)
iNKT = ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=iNKT, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
itm2 = ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=IntermediateMono_2, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
classical1 = ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=ClassicalMono_1, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
classical4 = ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=ClassicalMono_4, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
iNKToutlierrm = ggplot(PropMFIMetaData[PropMFIMetaData$iNKT<1.25,],aes(x=TimeToNextRelapse,y=iNKT, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
age = ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=Age, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
png(filename = 'Output/Plots/TimeToNextRelapseCorrelationPlotGridPerSex.png',
    width = 12,
    height = 8,
    units = 'in',
    res = 600)
ggarrange(plotlist = list(iNKT,itm2,classical1,iNKToutlierrm,classical4,age), common.legend = T)
dev.off()
png(filename = 'Output/Plots/TimeToNextRelapsePerTreatment.png',
    width = 12,
    height = 8,
    units = 'in',
    res = 600)
ggplot(PropMFIMetaData[!is.na(PropMFIMetaData$TreatmentSampling),], aes(x=TreatmentSampling, y=TimeToNextRelapse,fill=TreatmentSampling)) +  
  geom_violin(width=0.9,position = position_dodge(width = 0.8)) + 
  stat_summary(fun=median.quartile, position = position_dodge(width=0.8), geom = 'line') + 
  stat_summary(fun="median",size=1, geom="point",position = position_dodge(width=0.8)) +
  geom_quasirandom(size = 0.5,dodge.width = 0.8, varwidth = TRUE, alpha=0.3) + 
  stat_compare_means(size=7)+ylab('Time to next relapse (years)') + 
  xlab(NULL)+ theme_bw()+theme(text=element_text(size = 16), axis.text.x = element_text(size=14, angle = 45, hjust = 1), legend.text = element_text(size=14)) + 
  scale_fill_npg()
dev.off()

#Cor with spearman
ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=iNKT, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F,method = 'spearman')+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=IntermediateMono_1))+geom_point()+stat_cor(size=5, show.legend = F, method='spearman')+geom_smooth(method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=ClassicalMono_1, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=ClassicalMono_4, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()
ggplot(PropMFIMetaData[PropMFIMetaData$iNKT<1.25,],aes(x=TimeToNextRelapse,y=iNKT, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()


#Make ratio of ClassicalMono_1/ClassicalMono_4
PropMFIMetaData$ClassicalMono_1_ClassicalMono_4_Ratio=PropMFIMetaData$ClassicalMono_1/PropMFIMetaData$ClassicalMono_4
PropMFIMetaData <- PropMFIMetaData %>% relocate(ClassicalMono_1_ClassicalMono_4_Ratio, .after = Basophils)
PropMFIMetaData$ClassicalMono_4_ClassicalMono_1_Ratio=PropMFIMetaData$ClassicalMono_4/PropMFIMetaData$ClassicalMono_1
PropMFIMetaData <- PropMFIMetaData %>% relocate(ClassicalMono_4_ClassicalMono_1_Ratio, .after = ClassicalMono_1_ClassicalMono_4_Ratio)
ggplot(PropMFIMetaData,aes(x=TimeToNextRelapse,y=ClassicalMono_2+ClassicalMono_4, col=Sex))+geom_point()+stat_cor(aes(col=Sex),size=5, show.legend = F)+geom_smooth(aes(fill=Sex),method='lm')+theme_minimal() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) + scale_color_npg()


######################
#' 5 - Plots of correlation for time to SPMS
######################
cor_TimeToSPMS = find_cor(PropMFIMetaData, 'TimeToSPMS', c(1:170))
cor_TimeToSPMS$p.signif = makeStars(cor_TimeToSPMS$pvals)
png(filename = 'Output/Plots/TimeToSPMSSignificantParametersCorrelation.png',
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
ggplot(cor_TimeToSPMS[cor_TimeToSPMS$pvals<0.05,],aes(est,reorder(col_cluster,est)))+
  geom_col(width=.75,col='black',aes(fill=-log10(pvals)))+
  scale_fill_gradientn(colors=brewer.pal(name = 'YlOrRd', n = 9), name="-log10(p value)")+
  xlab("Correlation (r) with time to SPMS (years)")+ylab("")+
  geom_text(aes(label = p.signif,  hjust = ifelse(est < 0, 1.1, -0.1)), size=10, vjust = 0.8)+
  theme_bw()+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18), legend.text = element_text(size=14), legend.title=element_text(size=16))       
dev.off()

ggplot(PropMFIMetaData,aes(x=MFI.B2M_Bcells_2,y=TimeToSPMS))+geom_point()+geom_text(aes(label=rownames(PropMFIMetaData)))+stat_cor(size=7,show.legend = F)+geom_smooth(method='lm')+theme_minimal() +theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18)) + scale_color_npg()


ggplot(PropMFIMetaData[PropMFIMetaData$Cohort=='Outcome5yrs',], aes(x=NEDA5yrs, y=MFI.S100A8A9_NonClassicalMono_2,fill=NEDA5yrs)) +  
  geom_violin(width=0.9,position = position_dodge(width = 0.8)) + 
  stat_summary(fun=median.quartile, position = position_dodge(width=0.8), geom = 'line') + 
  stat_summary(fun="median",size=1, geom="point",position = position_dodge(width=0.8)) +
  geom_quasirandom(size = 0.5,dodge.width = 0.8, varwidth = TRUE, alpha=0.3) + 
  stat_compare_means(comparisons=list(c('NEDA',NA), c('EDA',NA)),size=7,method = 't.test')+ylab('MFI.S100A8A9_NonClassicalMono_2') + 
  xlab(NULL) +theme_bw()+theme(text=element_text(size = 16), axis.text.x = element_text(size=14, angle = 45, hjust = 1), legend.text = element_text(size=14)) + 
  scale_fill_npg()
######################
#' 6 - Plots of correlation for EDSS Sampling
######################
cor_EDSS_Sampling = find_cor(PropMFIMetaData[PropMFIMetaData$MSTypeSampling=='RRMS',], 'EDSS_Sampling', c(1:170))
cor_EDSS_Sampling$p.signif = makeStars(cor_EDSS_Sampling$pvals)
png(filename = 'Output/Plots/EDSS_SamplingRRMSSignificantParametersCorrelation.png',
    width = 13,
    height = 8,
    units = 'in',
    res = 600)
ggplot(cor_EDSS_Sampling[cor_EDSS_Sampling$pvals<0.05,],aes(est,reorder(col_cluster,est))) +
  geom_col(width=.75,col='black',aes(fill=-log10(pvals))) +
  scale_fill_gradientn(colors=brewer.pal(name = 'YlOrRd', n = 9), name="-log10(p value)") +
  xlab("Correlation (r) with EDSS at sampling")+ylab("") +
  geom_text(aes(label = p.signif,  hjust = ifelse(est < 0, 1.1, -0.1)), size=10, vjust = 0.8) + 
  theme_bw()+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18), legend.text = element_text(size=14), legend.title=element_text(size=16))
dev.off()
ggplot(PropMFIMetaData[PropMFIMetaData$MSTypeSampling=='RRMS',],aes(x=EDSS_Sampling,y=Age))+geom_point()+stat_cor(size=7,show.legend = F)+geom_smooth(method='lm')+theme_minimal() +theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18)) + scale_color_npg()
ggplot(PropMFIMetaData[PropMFIMetaData$MSTypeSampling=='RRMS',],aes(x=TimeToSPMS,y=TimeSinceOnset))+geom_point()+stat_cor(size=7,show.legend = F)+geom_smooth(method='lm')+theme_minimal() +theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18)) + scale_color_npg()
ggplot(PropMFIMetaData[PropMFIMetaData$MSTypeSampling=='RRMS',],aes(x=EDSS_Sampling,y=TimeToSPMS))+geom_point()+stat_cor(size=7,show.legend = F)+geom_smooth(method='lm')+theme_minimal() +theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18)) + scale_color_npg()

######################
#' 7 - Plots of correlation for delta EDSS
######################
cor_dEDSS = find_cor(PropMFIMetaData, 'dEDSS_5yrs', c(1:170))
cor_dEDSS$p.signif = makeStars(cor_dEDSS$pvals)
png(filename = 'Output/Plots/dEDSS_5yrsSignificantParametersCorrelation.png',
    width = 13,
    height = 8,
    units = 'in',
    res = 600)
ggplot(cor_dEDSS[cor_dEDSS$pvals<0.05,],aes(est,reorder(col_cluster,est))) +
  geom_col(width=.75,col='black',aes(fill=-log10(pvals))) +
  scale_fill_gradientn(colors=brewer.pal(name = 'YlOrRd', n = 9), name="-log10(p value)") +
  xlab("Correlation (r) with 5-year dEDSS")+ylab("") +
  geom_text(aes(label = p.signif,  hjust = ifelse(est < 0, 1.1, -0.1)), size=10, vjust = 0.8) + 
  theme_bw()+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 18), legend.text = element_text(size=14), legend.title=element_text(size=16))
dev.off()

