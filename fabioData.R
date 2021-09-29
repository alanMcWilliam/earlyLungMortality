
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(summarytools)


######################################################################


lungFabio <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/Lung_pathology.csv')
lungFabioCFS <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/clean_cfs.csv')
lungFabioStats <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/stats.csv')
lungFabioStatsExtra <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/extra_stats.csv')

lungFabioCheckExtra <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/seg_sanity_extra.csv')
lungFabioCheck <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/seg_sanity.csv')

lungFabioStats$ID <- str_remove(lungFabioStats$ID, '.npy')
lungFabioStats <- lungFabioStats %>%
  rename(patient_id = "ID",
         Area = "SM.Area",
         Density = "SM.Density")

lungFabioStatsExtra$ID <- str_remove(lungFabioStatsExtra$ID, '.npy')
lungFabioStatsExtra <- lungFabioStatsExtra %>%
  rename(patient_id = "ID",
         Area = "SM.Area",
         Density = "SM.Density")

lungFabioStats_combine <- rbind(lungFabioStats, lungFabioStatsExtra)
lungFabio <- merge(lungFabio, lungFabioStats_combine, by = 'patient_id')


View(lungFabio)
