# MetricEval_TestFunctions.r
# Purpose: Use this code to submit data to the functions in MetricEvaluation_General_forDistribution.r
# Created by Richard Mitchell 12/19/2022
####################################################################################################

library(plyr)
library(dplyr)
library(Hmisc)

setwd('C:/Users/rmitch11/OneDrive - Environmental Protection Agency (EPA)/VA MMI Analysis')

library('readxl')

#Reading in separate sheets from aquametBio_Final_03082022.xlsx file

VA_Site <- read_excel('aquametBio_final_03082022.xlsx', sheet='BenSamps')

Comp <- read_excel('aquametBio_final_03082022.xlsx', sheet='Composition')

Rich <- read_excel('aquametBio_final_03082022.xlsx', sheet='Richness')

Toler <- read_excel('aquametBio_final_03082022.xlsx', sheet='Tolerance')

Dominance <- read_excel('aquametBio_final_03082022.xlsx', sheet='Dominance')

Habit <- read_excel('aquametBio_final_03082022.xlsx', sheet='Habit')

FFG <- read_excel('aquametBio_final_03082022.xlsx', sheet='FFG')

library('tidyverse')

# Merge the various sheets

df_list <- list(VA_Site, Comp, Rich, Toler, Dominance, Habit, FFG)

df_list %>% reduce(full_join, by='BENSAMPID')

Merge <- merge(VA_Site, Comp, by='BENSAMPID')

Mer2 <- merge(Merge, Rich, by='BENSAMPID')

Mer3 <- merge(Mer2, Toler, by='BENSAMPID')

Mer4 <- merge(Mer3, Dominance, by='BENSAMPID')

Mer5 <- merge(Mer4, Habit, by='BENSAMPID')

Mer6 <- merge(Mer5, FFG, by='BENSAMPID')

Comb_Metrics <- Mer6

Comb_Metrics <- as.data.frame(Comb_Metrics)

Comb_Metrics[20:239][is.na(Comb_Metrics[20:239])] <- 0

write.csv(Comb_Metrics, 
          'C:/Users/rmitch11/OneDrive - Environmental Protection Agency (EPA)/VA MMI Analysis/VA_Comb_Met_New.csv'
          , row.names = F)

# I created separate seasonal files from the VA_Comb_Met_New.csv file; I did that in excel outside of R but you could also do in that in R as well.

dfIn <- read_excel('VA_MET_Coastal_Fall.xlsx', sheet='Sheet1') #Each season/bioregion dataframe needs to be run individually. 

dfIn.first <- dfIn

# I haven't run these yet, but this would be where I would create the Rep 1 sample Dataframe to run in the function.
# Need to create the individual season/bioregion files.
dfIn.Second <- read_excel('VA_MET_FALL_MOUNTAIN_V1_ONLY.xlsx', sheet='Sheet1') 

# The code below are the function calls for the various functions in Karen Blocksome's 'MettricEvaluation_General_forDistribution' code.
# Use only one visit per site. I didn't do this initially, so this will have to be rerun for the various seasonal/bioregion analsyes.
# Also, I did not run the S/N function becuase there are no sites that were sampled twice in the same season, but not on the sample date.
# We may want to consider running the S/N analysis on the multiple RepNum sites, but that will be a different analysis than what we did for NRSA and NLA.
test.rg <- rangeTest(dfIn.first,percVars=c('PTAX','PIND','%'),
                     idVars=c('BENSAMPID','StationID','RefStress','RepNum','Collection Date','Season','BioRegion'),
                     quantZero=0.75,
                     pass_=0.5,
                     quantRange=0.75,
                     quantRange.perc=15,
                     quantRange.oth=1/3)

test.red <- redTest(dfIn.first,c('BENSAMPID','StationID','RefStress','RepNum','Collection Date','Season','BioRegion'),0.8) # Skip this?

test.comp <- compTest(dfIn.first, idVars = c('BENSAMPID','StationID','RefStress','RepNum','Collection Date','Season','BioRegion'), refVar = 'RefStress','Ref','Str')

write.csv(test.comp,'C:/Users/rmitch11/OneDrive - Environmental Protection Agency (EPA)/VA MMI Analysis/K-W & BoxScore Comparison_Coastal_Fall.csv', row.names = F)

write.csv(test.rg,'C:/Users/rmitch11/OneDrive - Environmental Protection Agency (EPA)/VA MMI Analysis/Range_Test_Coastal_Fall.csv', row.names = F)

write.csv(test.red,'C:/Users/rmitch11/OneDrive - Environmental Protection Agency (EPA)/VA MMI Analysis/Redundancy_Analysis_Coastal_Fall.csv', row.names = F)