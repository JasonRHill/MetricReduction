# MetricEval_TestFunctions.r
# Purpose: Use this code to submit data to the functions in MetricEvaluation_General_forDistribution.r
# Created by Richard Mitchell 12/19/2022

###Edited by JRH 3/20/2023 in R 4.1.2
###note first loaded Karen Blocksom Functions before running this code:
###
# MetricEvaluation_Veg.R
## Created to run various screening tests on vegetation metrics and output results to Excel file
# 
# Author: kblockso
# 2/11/2014 Created by KAB
# 3/25/2014 Modified by KAB to add code to compare reference and trashed via boxplots and in Kruskal-Wallis test
# 4/14/2014 Modified by KAB to use updated REF_NWCA instead of RT_NWCA variable
# 10/15/2020 Updated by KAB to remove use of reshape2
##################################################################################################################
require(lme4)
require(Hmisc)
require(plyr)
require(dplyr)
require(gtools)
require(ggplot2)
require(tidyr)
library(tidyverse)

rangeTest <- function(dfIn,percVars,idVars,quantZero=0.75,pass_=NULL,quantRange=1,
                      quantRange.perc=15,quantRange.oth=1/3){
  # Purpose: To identify metrics with limited range or those that are highly skewed such that most of the values 
  #           are a single value. Allows user-specified parameters to perform these evaluations, with default 
  #           values provided.
  # Inputs:
  # dfIn - This data frame is in wide format with one row per sample and assumed to contain only numeric metrics. 
  #       It should also contain only the first visit to each site. If there are calibration and validation 
  #       subsets, this df should only contain calibration samples.
  #
  # percVars - A character vector containing either full names or partial names that clearly 
  #         identify the metrics that are percentages
  #
  # idVars - A character vector containing any variables that identify samples and are not metrics, 
  #         making it wise to drop any variables from the input data frame that are not necessary.
  #
  # quantZero - If this proportion of values is equal to the minimum value (typically 0), 
  #           metric fails range test. If the upper proportion of values is equal to the maximum value, 
  #           the metric also fails. Value should be a number between 0 and 1.
  #
  # pass_ - If the possibility of a partial pass is desired, provide an alternate value to quantZero which 
  #           is lower than quantZero. If this lower proportion of samples is equal to the minimum or the upper 
  #           proportion equal to this amount is the same as the max, a value of PASS- is assigned. If left blank, 
  #           this part is not performed. 
  #
  # quantRange - This value determines what upper proportion of sites is used to evaluate metric range. 
  #           For example, a value of 0.8 considers the range of values in the max - 20th percentile. 
  #           To examine the whole range, use a value of 1.
  #           Value should be a number between 0 and 1. 
  #
  # quantRange.perc - This numeric value represents a discrete value to which ranges of percentage metrics 
  #           are compared. For example, if quantRange=0.75 and quantRange.perc=15, if the range of the 
  #           upper 75% of values is less than 15, the metric fails.
  #
  # quantRange.oth - This numeric proportion of the full range is used as a point of comparison for the 
  #           range determined using quantRange. This is used for non-percentage metrics like richness. 
  #           For example, the default value of 1/3 multiplies the 
  #           full range by 1/3. For a quantRange of 0.75 and quantRange.perc of 1/3, 
  #           the difference between the 25th percentile and the max is compared to the max/3. 
  #           If the range is less than max/3, the metric fails. Value should be a number between 0 and 1.
  #
  #################################################################################################################
  metNames <- names(dfIn)[names(dfIn) %nin% idVars]
  percVars.1 <- paste(percVars,collapse='|')
  percNames <- grep(percVars.1,names(dfIn),value=TRUE)
  
  # Check ranges of values of input arguments
  if(quantZero<0|quantZero>1){
    return(print("quantZero value needs to be between 0 and 1"))    
  }
  if(quantRange<0|quantRange>1){
    return(print("quantRange value needs to be between 0 and 1"))
  }
  if(quantRange.oth<0|quantRange.oth>1){
    return(print("quantRange.oth value needs to be between 0 and 1"))
  }
  if(!is.null(pass_) & pass_<0|pass_>quantZero){
    return(print("pass_ value needs to be > 0 and < quantZero"))
  }
  
  # Now melt input data frame
  inLong <- tidyr::pivot_longer(dfIn, cols=names(dfIn)[names(dfIn) %nin% idVars], names_to='variable',
                                values_drop_na=TRUE) %>%
    dplyr::filter(!is.infinite(value)) %>% 
    mutate(value=as.numeric(value))
  
  # First get quantiles based on inputs
  rgQAll <- ddply(inLong,c('variable'),summarise,p0=min(value,na.rm=T),p100=max(value,na.rm=T)
                  ,plower=quantile(value,probs=(1-quantZero),na.rm=T),pupper=quantile(value,probs=(quantZero),na.rm=T),
                  prob.lower.rg=quantile(value, probs=(1-quantRange), na.rm=T))
  
  if(!is.null(pass_)){
    rgQpass_ <- ddply(inLong,c('variable'),summarise,pmid=quantile(value,probs=pass_,na.rm=T))
    rgQAll <- merge(rgQAll,rgQpass_,by='variable')
  }else{
    rgQAll <- mutate(pmid=NA)
  }
  
  # Range test - apply percentiles calculated above to test range of metrics
  rgTestAll <- mutate(rgQAll,zeroTest=ifelse(pupper==p0|plower==p100|pupper==0,'FAIL',
                                             ifelse(!is.null(pass_) & pmid==p0,'PASS-','PASS')),
                      rgLim=ifelse((variable %in% percNames & (p100-prob.lower.rg)<quantRange.perc)|
                                     (variable %nin% percNames & (p100-prob.lower.rg)<((p100-p0)*quantRange.oth)),'FAIL','PASS'))
  
  rgOutAll <- mutate(rgTestAll,RANGE_TEST=ifelse(zeroTest=='FAIL'|rgLim=='FAIL','FAIL'
                                                 ,ifelse(zeroTest=='PASS-','PASS-','PASS')),METRIC=as.character(variable)) %>%
    select(-variable)
  
  return(rgOutAll)
  
}



## Signal-to-noise Test
snTest <- function(dfIn,idVars.samp,idVars.site,year='YEAR'){
  # dfIn - This data frame is in wide format with one row per sample and assumed to contain only numeric metrics. 
  #       It should also contain all visits to each site, and at least a subset of sites must have multiple visits. 
  #       If there are calibration and validation subsets, this df should only contain calibration samples.
  #       Only identifying variables and metrics should be included in this data frame, and only numeric metrics 
  #       can be included.
  # 
  # idVars.samp - a character vector containing variables that identify individual samples. 
  #
  # idVars.site - a string containing variable name that identifies sites. This cannot be the same as or a subset 
  #       of variables in idVars.samp.
  #
  # year - string containing name of Year variable if sites are revisited across years (as well as within year), 
  #   default is 'YEAR'. Set to NULL if no samples across years.
  #############################################################################################################################
  options(warn=2)
  
  # Do some error checking first
  if(idVars.site %in% idVars.samp){
    return(print("idVars.site CANNOT be a part of idVars.samp"))
  }
  # Make sure there are some repeated sites
  if(nrow(dfIn[duplicated(dfIn[,idVars.site]),])<1){
    return(print("You need multiple visits for at least 1 site"))
  }  
  print(c('Number of revisits:',nrow(dfIn[duplicated(dfIn[,idVars.site]),])))
  if(!is.null(year)){
    inLong <- tidyr::pivot_longer(dfIn, cols=names(dfIn)[names(dfIn) %nin% c(idVars.samp, idVars.site, year)],
                                  names_to='variable', values_drop_na=TRUE) %>%
      dplyr::filter(!is.infinite(value)) %>%
      mutate(variable=as.character(variable),value=as.numeric(value)) 
    names(inLong)[names(inLong)==idVars.site] <- 'site'
    names(inLong)[names(inLong)==year] <- 'year'
    
  }else{
    inLong <- tidyr::pivot_longer(dfIn, cols=names(dfIn)[names(dfIn) %nin% c(idVars.samp, idVars.site)],
                                  names_to='variable', values_drop_na=TRUE) %>%
      dplyr::filter(!is.infinite(value)) %>%
      mutate(variable=as.character(variable),value=as.numeric(value)) 
    names(inLong)[names(inLong)==idVars.site] <- 'site' 
    
  }
  
  # create vector of metric names
  parList <- unique(inLong$variable)
  
  ## Signal-to-noise ratio
  # Create empty data frame to accept output for each metric
  snOut <- data.frame(METRIC=character(),SIGNAL=numeric(),NOISE=numeric(),SN_RATIO=numeric(),COM=character(),stringsAsFactors=FALSE)
  
  # For each metric in parList, run a linear mixed-effects model with SITE_ID as a random effect
  for(i in 1:length(parList)){
    inMet <- subset(inLong,variable==parList[i])
    # print(parList[i])					
    # Run model
    if(!is.null(year)){
      sn <- try(lmer(value~year + (1|year:site),inMet),
                silent=TRUE)
      
      if(class(sn)=='try-error'){
        sn <- try(lmer(value~year + (1|year:site),inMet, control= lmerControl(optimizer = "bobyqa", 
                                                                              optCtrl = list(maxfun = 100000))),
                  silent=TRUE)
      }
    }else{
      sn <- try(lmer(value~(1|site),inMet),silent=TRUE)
    }
    
    # If model output is error, send to output data frame
    if(class(sn)=='try-error'){
      StoN <- data.frame(METRIC=parList[i],SIGNAL=NA,NOISE=NA,SN_RATIO=NA,COM=sn[1],stringsAsFactors=FALSE)
    }else{
      # If model output value, determine signal and noise by extracting variance components due to SITE_ID and error
      varcomp <- VarCorr(sn)
      if(!is.null(year)){
        StoN <- data.frame(METRIC=parList[i],SIGNAL=round(varcomp$'year:site'[1],2), 
                           NOISE=round(attr(varcomp,"sc")^2,2)
                           ,SN_RATIO=round(varcomp$'year:site'[1]/(attr(varcomp,"sc")^2),2),COM=NA,
                           stringsAsFactors=FALSE)
      }else{
        StoN <- data.frame(METRIC=parList[i], SIGNAL= round(varcomp$site[1],2), NOISE=round(attr(varcomp,"sc")^2,2)
                           ,SN_RATIO=round(varcomp$site[1]/(attr(varcomp,"sc")^2),2),COM=NA,stringsAsFactors=FALSE)
      }
    }			
    snOut <- rbind(snOut,StoN)
  }
  return(snOut)
  
}



redTest <- function(dfIn,idVars,cutoff){
  # dfIn - This is a wide data frame with one row per sample. It is assumed to only contain one sample per site,
  #         so if there are revisits, keep only the first visit. It is assumed that all metrics are numeric, so 
  #         remove any that are not.
  # idVars - This is a character vector with any non-metric variables in the data frame
  # cutoff - Minimum absolute correlation value at which metrics are considered redundant
  #####################################################################################################################
  
  # Melt input dataset by all variables in data frame that are not numeric metrics (or remove those variables)
  corrIn <- tidyr::pivot_longer(dfIn, cols = names(dfIn)[names(dfIn) %nin% idVars], 
                                names_to='variable', values_drop_na=TRUE) %>%
    # reshape2::melt(dfIn,id.vars=idVars,na.rm=T) %>% 
    filter(!is.infinite(value)) %>%
    mutate(value=as.numeric(value),variable=as.character(variable))
  # Create empty data frame 
  corrOut <- data.frame(METRIC=character(),REDUND_MET=character(),stringsAsFactors=FALSE)
  
  metList <-unique(corrIn$variable)
  for(i in 1:length(metList)){
    print(i)
    # Make sure to alter to exclude the id variables used in melting the input data
    corMet <- cor(subset(corrIn,variable==metList[i],select='value')
                  ,subset(dfIn,select=names(dfIn) %nin% c(idVars,metList[i]))
                  ,method="pearson")
    redMet <- data.frame(METRIC=attr(corMet,"dimnames")[[2]],R=corMet[1:length(corMet)]
                         ,stringsAsFactors=FALSE) %>% 
      subset(abs(R)>=cutoff,select='METRIC')
    redList <- data.frame(METRIC=metList[i],REDUND_MET=paste(redMet$METRIC,collapse=","),stringsAsFactors=FALSE)
    corrOut <- rbind(corrOut,redList)
  }
  return(corrOut)
}



compTest <- function(dfIn,idVars,refVar,least,most){
  # This function compares the most and least disturbed using a Kruskal-Wallis test because metrics are typically somewhat
  #       skewed. It also performs a box plot comparison using scoring developed by Tetra Tech for Florida SCI.
  #
  # dfIn - This is a wide data frame with one row per sample. It is assumed to only contain one sample per site,
  #         so if there are revisits, keep only the first visit. It is assumed that all metrics are numeric, so 
  #         remove any that are not. This data frame should include variables to identify samples and a single variable 
  #         for level of disturbance. The data frame should include only the most and least disturbed classes. Some 
  #         consolidation of classes may be desirable, depending on the situation.
  #
  # idVars - This is a character vector with variables to identify samples in the data frame. 
  #
  # refVar - string identifying the variable with disturbance condition for each site.
  #
  # least - string representing Value of refVar indicating least disturbed condition
  #
  # most - string representing value of refVar indicating most disturbed condition
  ###########################################################################################################################
  inLong <- tidyr::pivot_longer(dfIn, cols = names(dfIn)[names(dfIn) %nin% c(idVars, refVar)], 
                                names_to='variable', values_drop_na=TRUE) %>%
    subset(eval(as.name(refVar)) %in% c(least,most)) %>% 
    mutate(variable=as.character(variable),value=as.numeric(value))
  
  # Create empty data frame to accept test output
  compOut <- data.frame(METRIC=character(),KW_stat=numeric(),KW_pval=numeric(),BoxScore=integer(),stringsAsFactors=FALSE)
  # Create list of parameters
  metList <- as.character(unique(inLong$variable))
  
  # For each metric, first run a Kruskal-Wallis test comparing best and worst categories
  # Then determine quantiles to simulate boxplot comparisons as used by Tetra Tech
  for(i in 1:length(metList)){
    # print(metList[i])
    ## Obtain K-W p-value
    kw <- kruskal.test(value~factor(eval(as.name(refVar))),data=subset(inLong,variable==metList[i]))
    
    ## Use interquartile ranges to simulate comparison of boxplots and use scoring from Tetra Tech
    quants <- ddply(subset(inLong,variable==metList[i]),c(refVar),summarise,p25=quantile(value,probs=0.25)
                    ,p75=quantile(value,probs=0.75),median=quantile(value,probs=0.5))
    
    quants.1 <- data.frame(METRIC=metList[i],
                           tidyr::pivot_longer(quants, cols=names(quants)[names(quants) %nin% refVar], names_to='variable')) %>%
      #reshape2::melt(quants,id.vars=c(refVar))) %>% 
      mutate(variable = paste(eval(as.name(refVar)), variable, sep='_')) %>%
      tidyr::pivot_wider(id_cols='METRIC', names_from='variable')
    #dcast(METRIC~eval(as.name(refVar))+variable)
    
    quants.2 <- mutate(quants.1,munder=ifelse(Str_median<Ref_p25,1,0),mover=ifelse(Str_median>Ref_p75,1,0)
                       ,lunder=ifelse(Ref_median<Str_p25,1,0)
                       ,lover=ifelse(Ref_median>Str_p75,1,0)
                       ,overlap=ifelse((munder==1 & lover==1 & Str_p75<Ref_p25)|(mover==1 & lunder==1 & Ref_p75<Str_p25),1,0)
                       ,BoxScore=sum(munder,mover,lunder,lover,overlap))
    
    tempComp <- data.frame(METRIC=metList[i],KW_stat=round(kw$statistic,2),KW_pval=round(kw$p.value,4),BoxScore=quants.2$BoxScore
                           ,stringsAsFactors=FALSE)
    compOut <- rbind(compOut,tempComp)
  }
  return(compOut)
}

###Start With All Data
dfIn <- read.csv("Data/VA_Comb_Met_New_All.csv")


###lou wants to know all the 1's so we can remove from the data set to start with

test.red.all.lou <- redTest(dfIn,c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),1)

write.csv(test.red.all.lou,'Output/Redundancy_Analysis_All_Lou.csv', row.names = F)

# The code below are the function calls for the various functions in Karen Blocksome's 'MettricEvaluation_General_forDistribution' code.
# This first run includes all metrics for all sites in all regions. Super high level, will goto non-coastal area all, then non-coastal season.

##Then will look at Richards will bind and make some recommendations.


# Also, I did not run the S/N function because there are no sites that were sampled twice in the same season, but not on the sample date.
# We may want to consider running the S/N analysis on the multiple RepNum sites, but that will be a different analysis than what we did for NRSA and NLA.

test.rg.all <- rangeTest(dfIn,percVars=c('PTAX','PIND','%'),
                     idVars=c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),
                     quantZero=0.75,
                     pass_=0.5,
                     quantRange=0.75,
                     quantRange.perc=15,
                     quantRange.oth=1/3)

test.red.all <- redTest(dfIn,c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),0.8)

test.comp.all <- compTest(dfIn, idVars = c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'), refVar = 'RefStress','Ref','Str')


write.csv(test.rg.all,'Output/Range_Test_All.csv', row.names = F)

write.csv(test.red.all,'Output/Redundancy_Analysis_All.csv', row.names = F)

write.csv(test.comp.all,'Output/K-W_BoxScore_Comparison_All.csv', row.names = F)


###Start With All Data in Non-Coastal Areas
dfIn1 <- read.csv("Data/VA_Comb_Met_New_All_NoCoast.csv")


# The code below are the function calls for the various functions in Karen Blocksome's 'MettricEvaluation_General_forDistribution' code.
# This first run includes all metrics for all sites in all regions. Super high level, will goto non-coastal area all, then non-coastal season.

##Then will look at Richards will bind and make some recommendations.


# Also, I did not run the S/N function because there are no sites that were sampled twice in the same season, but not on the sample date.
# We may want to consider running the S/N analysis on the multiple RepNum sites, but that will be a different analysis than what we did for NRSA and NLA.

test.rg.all.nocoast <- rangeTest(dfIn1,percVars=c('PTAX','PIND','%'),
                         idVars=c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),
                         quantZero=0.75,
                         pass_=0.5,
                         quantRange=0.75,
                         quantRange.perc=15,
                         quantRange.oth=1/3)

test.red.all.nocoast <- redTest(dfIn1,c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),0.8)

test.comp.all.nocoast <- compTest(dfIn1, idVars = c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'), refVar = 'RefStress','Ref','Str')


write.csv(test.rg.all.nocoast,'Output/Range_Test_All_NoCoast.csv', row.names = F)

write.csv(test.red.all.nocoast,'Output/Redundancy_Analysis_All_NoCoast.csv', row.names = F)

write.csv(test.comp.all.nocoast,'Output/K-W_BoxScore_Comparison_AllNoCoast.csv', row.names = F)


###Start With All Data in Just Coastal Areas
dfIn2 <- read.csv("Data/VA_Comb_Met_New_All_JustCoast.csv")


# The code below are the function calls for the various functions in Karen Blocksome's 'MettricEvaluation_General_forDistribution' code.
# This first run includes all metrics for all sites in all regions. Super high level, will goto non-coastal area all, then non-coastal season.

##Then will look at Richards will bind and make some recommendations.


# Also, I did not run the S/N function because there are no sites that were sampled twice in the same season, but not on the sample date.
# We may want to consider running the S/N analysis on the multiple RepNum sites, but that will be a different analysis than what we did for NRSA and NLA.

test.rg.all.justcoast <- rangeTest(dfIn2,percVars=c('PTAX','PIND','%'),
                                 idVars=c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),
                                 quantZero=0.75,
                                 pass_=0.5,
                                 quantRange=0.75,
                                 quantRange.perc=15,
                                 quantRange.oth=1/3)

test.red.all.justcoast <- redTest(dfIn2,c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),0.8)

test.comp.all.justcoast <- compTest(dfIn2, idVars = c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'), refVar = 'RefStress','Ref','Str')


write.csv(test.rg.all.justcoast,'Output/Range_Test_All_JustCoast.csv', row.names = F)

write.csv(test.red.all.justcoast,'Output/Redundancy_Analysis_All_JustCoast.csv', row.names = F)

write.csv(test.comp.all.justcoast,'Output/K-W_BoxScore_Comparison_JustCoast.csv', row.names = F)


###Start With All Data in Non_Coastal Areas Spring
dfIn3 <- read.csv("Data/VA_Comb_Met_New_All_NoCoastSpring.csv")


# The code below are the function calls for the various functions in Karen Blocksome's 'MettricEvaluation_General_forDistribution' code.
# This first run includes all metrics for all sites in all regions. Super high level, will goto non-coastal area all, then non-coastal season.

##Then will look at Richards will bind and make some recommendations.


# Also, I did not run the S/N function because there are no sites that were sampled twice in the same season, but not on the sample date.
# We may want to consider running the S/N analysis on the multiple RepNum sites, but that will be a different analysis than what we did for NRSA and NLA.

test.rg.all.nocoast.spring <- rangeTest(dfIn3,percVars=c('PTAX','PIND','%'),
                                   idVars=c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),
                                   quantZero=0.75,
                                   pass_=0.5,
                                   quantRange=0.75,
                                   quantRange.perc=15,
                                   quantRange.oth=1/3)

test.red.all.nocoast.spring <- redTest(dfIn3,c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),0.8)

test.comp.all.nocoast.spring <- compTest(dfIn3, idVars = c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'), refVar = 'RefStress','Ref','Str')


write.csv(test.rg.all.nocoast.spring,'Output/Range_Test_All_NoCoast_Spring.csv', row.names = F)

write.csv(test.red.all.nocoast.spring,'Output/Redundancy_Analysis_All_NoCoast_Spring.csv', row.names = F)

write.csv(test.comp.all.nocoast.spring,'Output/K-W_BoxScore_Comparison_NoCoast_Spring.csv', row.names = F)


###Start With All Data in Non_Coastal Areas Fall
dfIn4 <- read.csv("Data/VA_Comb_Met_New_All_NoCoastFall.csv")


# The code below are the function calls for the various functions in Karen Blocksome's 'MettricEvaluation_General_forDistribution' code.
# This first run includes all metrics for all sites in all regions. Super high level, will goto non-coastal area all, then non-coastal season.

##Then will look at Richards will bind and make some recommendations.


# Also, I did not run the S/N function because there are no sites that were sampled twice in the same season, but not on the sample date.
# We may want to consider running the S/N analysis on the multiple RepNum sites, but that will be a different analysis than what we did for NRSA and NLA.

test.rg.all.nocoast.fall <- rangeTest(dfIn4,percVars=c('PTAX','PIND','%'),
                                        idVars=c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),
                                        quantZero=0.75,
                                        pass_=0.5,
                                        quantRange=0.75,
                                        quantRange.perc=15,
                                        quantRange.oth=1/3)

test.red.all.nocoast.fall <- redTest(dfIn4,c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),0.8)

test.comp.all.nocoast.fall <- compTest(dfIn4, idVars = c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'), refVar = 'RefStress','Ref','Str')


write.csv(test.rg.all.nocoast.fall,'Output/Range_Test_All_NoCoast_Fall.csv', row.names = F)

write.csv(test.red.all.nocoast.fall,'Output/Redundancy_Analysis_All_NoCoast_Fall.csv', row.names = F)

write.csv(test.comp.all.nocoast.fall,'Output/K-W_BoxScore_Comparison_NoCoast_Fall.csv', row.names = F)




###Check Greg Collection Concern
dfIn5 <- read.csv("Data/VA_Comb_Met_New_All_MACS.csv")


# The code below are the function calls for the various functions in Karen Blocksome's 'MettricEvaluation_General_forDistribution' code.
# This first run includes all metrics for all sites in all regions. Super high level, will goto non-coastal area all, then non-coastal season.

##Then will look at Richards will bind and make some recommendations.


# Also, I did not run the S/N function because there are no sites that were sampled twice in the same season, but not on the sample date.
# We may want to consider running the S/N analysis on the multiple RepNum sites, but that will be a different analysis than what we did for NRSA and NLA.

test.rg.all.MACS <- rangeTest(dfIn5,percVars=c('PTAX','PIND','%'),
                                      idVars=c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),
                                      quantZero=0.75,
                                      pass_=0.5,
                                      quantRange=0.75,
                                      quantRange.perc=15,
                                      quantRange.oth=1/3)

test.red.all.MACS <- redTest(dfIn5,c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'),0.8)

test.comp.all.MACS <- compTest(dfIn5, idVars = c('BENSAMPID','StationID','RefStress','RepNum','CollectionDate','Gradient','Season','BioRegion'), refVar = 'RefStress','Ref','Str')


write.csv(test.rg.all.MACS,'Output/Range_Test_All_MACS.csv', row.names = F)

write.csv(test.red.all.MACS,'Output/Redundancy_Analysis_All_MACS.csv', row.names = F)

write.csv(test.comp.all.MACS,'Output/K-W_BoxScore_Comparison_MACS.csv', row.names = F)




