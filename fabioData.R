
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr) ## 
library(summarytools)
library(caret)

#install.packages('caret')

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


lungFabio <- lungFabio %>%
  filter(!patient_id %in% lungFabioCheck$ID) %>%
  filter(!patient_id %in% lungFabioCheckExtra$ID)

View(lungFabio)


### select variables for models
#### remove patients treated prior to 2018 to avoid any potential overlap in data, but very unlikely

summary(lungFabio$Overall.Survival.Time..Days.) 

lungFabio$nintyDay <- lungFabio$Overall.Survival.Time..Days. <= 115
lungFabio$hundredEighty <- lungFabio$Overall.Survival.Time..Days. <= 200

summary(lungFabio$nintyDay)
summary(lungFabio$hundredEighty)

### need to consider what to do with the chemo data....;
lungFabioClean <- lungFabio %>%
  filter(as.Date(Mosaiq.Date.of.First.Fraction.Primary.Course) >= "2017-01-01") %>%
  filter(Primary.Intent == 'Radical') %>%
  filter(M.Stage == '0') %>%
  rename(Age = "Age.start.treatment",
         T.stageClean = "T.Stage",
         N.stageClean = "N.Stage",
         performance.status = "ecog_ps",
         gender = "sex",
         fullArea = "Area",
         halfDensity = "Density") %>%
  select(performance.status, gender, smoking, Pack.Years, comorbidities, Primary.Prescribed.Dose, Primary.Prescribed.Fractions, Primary.Intent, Overall.Survival.Time..Days., Death.Status, Age, T.stageClean, N.stageClean, M.Stage, fullArea, halfDensity, nintyDay, hundredEighty)

### convert to cm3 assuming sampled at ~1mm
lungFabioClean$fullAreaCC <- lungFabioClean$fullArea*0.01

lungFabioClean$T.stageClean <- as.factor(lungFabioClean$T.stageClean)
lungFabioClean$N.stageClean <- as.factor(lungFabioClean$N.stageClean)


stview(dfSummary(lungFabioClean))


#########################################################
### modeling building

lungFabioNew <- lungFabioClean %>%
  select(gender, Age, T.stageClean, N.stageClean, fullAreaCC, nintyDay)
lungFabioNew <- lungFabioNew[complete.cases(lungFabioNew),]

lungFabioNew$nintyDay <- as.numeric(lungFabioNew$nintyDay)
#lungFabioNew$nintyDay <- as.factor(lungFabioNew$nintyDay)


sarcNintyFabio <- glm(nintyDay~fullAreaCC + Age + T.stageClean + N.stageClean + gender, data = lungFabioClean, family=binomial(link='logit'))
summary(sarcNintyFabio)
summary(sarcNinty)

pred <- predict(sarcNinty, newdata = lungFabioNew)
#lungFabioNew$pred_glm = ifelse(pred > 0.5, "1", "0")
#lungFabioNew$pred_glm = as.factor(lungFabioNew$pred_glm)

#confusionMatrix(lungFabioNew$nintyDay, lungFabioNew$pred_glm)
#summary(pred)

table(pred)
table(lungFabioNew$nintyDay)
confusionMatrix(as.factor(pred), as.factor(lungFabioNew$nintyDay))


probabilities <- sarcNinty %>% predict(lungFabioNew, type = "response")
head(probabilities)

#contrasts(lungFabioNew$nintyDay)
predicted.classes <- ifelse(probabilities > 0.2, "1", "0")
head(predicted.classes)

mean(predicted.classes == lungFabioNew$nintyDay)
confusionMatrix(as.factor(predicted.classes), as.factor(lungFabioNew$nintyDay))



##### combine Fabio and standard data


sarcDensityFabio <- glm(lungFabioClean$nintyDay~lungFabioClean$Density)
summary(sarcDensityFabio)

sarcDensityFabio <- glm(lungFabioClean$hundredEighty~lungFabioClean$Area)
summary(sarcDensityFabio)

sarcDensityFabio2 <- glm(hundredEighty~Area + Age.start.treatment + sex + factor(T.Stage) + factor(N.Stage), data = lungFabioClean)
summary(sarcDensityFabio2)

