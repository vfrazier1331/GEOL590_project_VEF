install.packages("readxl")
library(readxl)
library(tidyverse)

anionstd_10 <- read_xls("data/ASII_10.xls")
print(anionstd_10)

head(anionstd_10)

results <- anionstd_10[42:48, 2:8] #pulls out results ... need to label columns

results_df <- data.frame(results) #creates a data frame out of the data we pulled out
print(results_df)
#renaming columns to be informative
names(results_df) <- c("anion", "retention_time.min", "area.uS*min", "height.uS", "relative_area.%", "relative_height.%", "concentration.mg/L")

run_info1.0 <- anionstd_10[3:9, 6:7] #pulls out information about the run
names(run_info1) <- c("parameter", "data") #renames columns
print(run_info1)

run_info2 <- anionstd_10[3:9, 1:3] #second set of info
run_info2.0 <- run_info2[,-2] #removes the empty column in second set
names(run_info2.0) <- c("parameter", "data")
print(run_info2.0)

run_info_prelim <- data.frame(run_info1.0) #makes run_info1.0 into dataframe. will add rows for run_info2.0
run_info_df <- rbind(run_info2.0, run_info1.0)
print(run_info_df)
glimpse(run_info_df)
str(run_info_df)
#we need to keep info with the data, so I'm going to create a list with both dataframes
output_list <- list(results_df, run_info_df)

