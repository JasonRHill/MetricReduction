
# MMI METRICS FORMATTING #
# This code serves to format the allMetrics dataframe to contain only those
# columns which are specific to Coast and NonCoast sampling.
# Author: Kansas Keeton
# Date Last Edited: 1/6/2024

library(tidyverse)
library(openxlsx)

# Create function to read all Excel sheets
readExcel <- function(sheetIndex){
  d <- read.xlsx("Final_Metrics_Jan2024_Jason_to_Kansas.xlsx", sheet = sheetIndex)
}
# create vector of sheet indices as input for function
sheetIndices <- c(1, 2, 3)
# map vector to function
file <- map(sheetIndices, readExcel)
# name list elements
names(file) <- c("allMetrics", "nonCoast", "Coast")
# list data frames to the global environment
list2env(file, globalenv())
# Create unique "Metric_Type" columns since R requires uniquely named cols
metricColumns <- names(allMetrics) %>%
  as_tibble() %>%
  mutate(newNames = if_else(value == "Metric_Type",
                            paste(lag(value), value, sep = "_"),
                            value))
# Preserve raw data by creating new dataframe
newMetrics <- allMetrics
# assign column names
names(newMetrics) <- metricColumns$newNames
# change column formats to match column names because Excel reads columns weird
names(newMetrics) <- gsub(".-.", " - ", names(newMetrics), fixed = TRUE)
names(newMetrics) <- gsub(".", " ", names(newMetrics), fixed = TRUE)
names(newMetrics) <- gsub("T0", "To", names(newMetrics), fixed = TRUE)
names(newMetrics) <- gsub("4 5", "4.5", names(newMetrics), fixed = TRUE)
names(newMetrics) <- gsub("6 5", "6.5", names(newMetrics), fixed = TRUE)
names(newMetrics) <- gsub("TN TP", "TN.TP", names(newMetrics), fixed = TRUE)
# Add in Coast Metric Types
Coast$Metric_Type <- paste(Coast$METRIC, "Metric_Type", sep = "_")
# Combine the column names
coastAllNames <- as.vector(t(Coast))
# select dataframe by metric names
coastMetrics <- newMetrics %>%
  select(DataSource:`Entered By`, all_of(coastAllNames)) %>% # select ID columns and Coast Metrics
  mutate(RepNum = ifelse(RepNum==0, 1, RepNum),
         Gradient = ifelse(is.na(Gradient), "Riffle", Gradient)) %>%
  filter(RepNum==1) %>% # remove all reps except 1
  filter(BioRegion=="Coast") %>%
  filter(!str_detect(BenSampID, regex("\\QAQC"))) # remove QAQC samples

# Add in nonCoast Metric Types
nonCoast$Metric_Type <- paste(nonCoast$MetricNonCoast, "Metric_Type", sep = "_")
# combine the column names
nonCoastAllNames <- as.vector(t(nonCoast))
# select dataframe by metric names
nonCoastMetrics <- newMetrics %>%
  select(DataSource:`Entered By`, all_of(nonCoastAllNames)) %>% # select ID columns and Coast Metrics
  mutate(RepNum = ifelse(RepNum==0, 1, RepNum),
         Gradient = ifelse(is.na(Gradient), "Riffle", Gradient)) %>%
  filter(RepNum==1) %>% # remove all reps except 1
  filter(BioRegion!="Coast") %>%
  filter(!str_detect(BenSampID, regex("\\QAQC"))) # remove QAQC samples

# Write to Excel file
list_of_metrics <- list("Raw" = allMetrics, "Coast_Metrics" = coastMetrics,
                        "NonCoast_Metrics" = nonCoastMetrics)
write.xlsx(list_of_metrics, file = paste0("Final_Metrics_Formatted_", Sys.Date(), ".xlsx"))
