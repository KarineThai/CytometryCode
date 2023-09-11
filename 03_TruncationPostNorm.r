########################################
#' Cut outlier events in normalized FCS - min and max gates done in FlowJo and saved in excel table
#' Generates files that are transformed and normalized for batch effects, exported as FCS and as RData flowset
#' Karine Thai - 2023/07/12, Updated on 2023/07/24
########################################

######################
#' 0 - load required packages
######################
library(ggplot2)
library(ggridges)
library(flowCore)
library(reshape2)
library(gridExtra)
library(umap)
library(RColorBrewer)
library(viridis)
library(flowCore)
library(readxl)
library(stringr)

######################
#' 1 - Read flowset and rename files that had wrong label, if needed
######################
fs_Norm = read.flowSet(path = 'Output/Normalized_FCS/', pattern="fcs$", truncate_max_range=T)
fs_Norm@frames$Norm_batch4_mspatient_stain_gb_mg_128_20140812.fcs@description$`$FIL`="batch4_mspatient_stain_gb_mg_202_20140812.fcs"
fs_Norm@frames$Norm_batch4_mspatient_stain_jl_mg_202_20140508.fcs@description$`$FIL`="Norm_batch4_mspatient_stain_jl_mg_185_20140508.fcs"
fs_Norm@frames$Norm_batch4_hlctrl_stain_cd_159_20210503.fcs@description$`$FIL`="Norm_batch4_hlctrl_stain_cp_159_20210503.fcs"
save(fs_Norm, file = 'Output/RData/fs_Norm.RData')

######################
#' 2 - Create a list of tables of expression where each list element is a sample
######################

load('Output/RData/fs_Norm.RData')
source("Functions/Create_list_expression_table_from_flowset.R")

regex1="(mspatient|hlctrl)"
regex2="([a-z]{2,3})_([0-9]{3,4})"
str_match(fs_Norm@phenoData@data$name, "(mspatient|hlctrl)")
Tables_exp=Create_list_expression_table_from_flowset(fs_Norm,regex1=regex1, regex2=regex2)
names(Tables_exp)
head(Tables_exp[[1]])
save(Tables_exp, file='Output/RData/Tables_exp.RData')

######################
#' 3 - Read table containing ranges of expression values for each marker
######################

Table_ParameterRange <- as.data.frame(read_excel("Input_Tables/Table_ParameterRange.xlsx",col_names = T))

######################
#' 4 - Create a new list of tables without events outside of the range
######################

Tables_exp_cut <- list()

# Iterate over each dataframe in 'Tables_exp'
for (i in seq_along(Tables_exp)) {
  # Get the name of the current dataframe
  df_name <- names(Tables_exp)[i]
  
  # Subset the corresponding rows based on the desired ranges
  filtered_df <- Tables_exp[[i]][apply(Tables_exp[[i]], 1, function(row) {
    # Check if each observation type falls within the desired range
    all(row[Table_ParameterRange$Param] >= Table_ParameterRange$Min &
          row[Table_ParameterRange$Param] <= Table_ParameterRange$Max)
  }), ]
  
  # Assign the filtered dataframe to the new list with the same dataframe name
  Tables_exp_cut[[df_name]] <- filtered_df
}

# Return the new list of filtered dataframes
Tables_exp_cut

######################
#' 5 - Verify the percentage of cells removed in each sample
######################

source('Functions/NumberCellsPerSampleTruncated.R')
SummaryTruncation=NumberCellsPerSampleTruncated(Tables_exp,Tables_exp_cut)
SummaryTruncation

######################
#' 6 - Exclude samples with more than 4% cells removed and save 
######################

TruncationOver4Percent=SummaryTruncation[SummaryTruncation$Percent_cut>4,]
ListTablesExpCut=Tables_exp_cut[names(Tables_exp_cut) %in% rownames(TruncationOver4Percent) == FALSE]
######################
#' 7 - Remove CD45 parameter because not needed for further analyses
######################

ListTablesExpCut <- lapply(ListTablesExpCut, function(df) df[, !names(df) %in% "CD45", drop = FALSE])
save(ListTablesExpCut, file = 'Output/RData/ListTablesExpCut.RData')

