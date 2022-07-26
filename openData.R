
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(summarytools)
library(tidyr)
library(survival)


######################################################################


lungOpen <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/openData/data/data/statistics/T12_statistics.csv')
lungOpenClinical <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/openData/NSCLCradiomicsClin.csv')

lungOpen$filename <- str_remove(lungOpen$filename, '.nii')
lungOpen <- lungOpen %>%
  rename(PatientID = "filename") %>%
  filter(sm_density > -20)

lungOpen <- merge(lungOpenClinical, lungOpen, by = 'PatientID')
lungOpen <- lungOpen%>%
  filter(clinical.T.Stage < 5)

lungOpen$sm_area <- lungOpen$sm_area*0.001
lungOpen$nintyDay <- lungOpen$Survival.time < 120

 
lungOpen$densityMedian <- lungOpen$sm_density < median(lungOpen$sm_density)
### 12.5 median from Manchester cohort

lungOpen$densityMedian <- lungOpen$sm_density < median(lungOpen$sm_density)


# + Age + T.stageClean + N.stageClean + gender
openDensity <- glm(nintyDay~ sm_density + age + gender + factor(clinical.T.Stage) + factor(Clinical.N.Stage), data = lungOpen, family="binomial")
summary(openDensity)

summary(lungOpen$sm_density)
summary(lungOpen$sm_area)
summary(lungOpen$nintyDay)

summary(lungOpen$age)

tapply(lungOpen$sm_density, lungOpen$gender, summary)
tapply(lungOpen$sm_area, lungOpen$gender, summary)


ggplot(data=lungOpen, aes(x = sm_density, color = gender)) + 
  geom_histogram(breaks=seq(-30,60, by = 2)) +
  labs(title = "", x = "Muscle density" ) +
  theme(panel.background = element_blank())

lungOpenMale <- lungOpen %>%
  filter(gender == "male")
lungOpenFemale <- lungOpen %>%
  filter(gender == "female")


summary(lungOpenMale$sm_density)
summary(lungOpenFemale$sm_density)


# median in Manchester M 12.8, F 12.2
# medin areas in Manchester patients M 51.5, F 35
lungOpenMale$densityMedianSS <- lungOpenMale$sm_density < median(lungOpenMale$sm_density)
lungOpenMale$densityMedianSS_Man <- lungOpenMale$sm_density < 15.7

lungOpenMale$densityAreaSS <- lungOpenMale$sm_area < median(lungOpenMale$sm_area)

lungOpenFemale$densityMedianSS <- lungOpenFemale$sm_density < median(lungOpenFemale$sm_density)
lungOpenFemale$densityMedianSS_Man <- lungOpenFemale$sm_density < 10

lungOpenFemale$densityAreaSS <- lungOpenFemale$sm_area < median(lungOpenFemale$sm_area)

lungOpenNew <- rbind(lungOpenMale, lungOpenFemale)

openDensityM <- glm(nintyDay~ sm_area + age + factor(clinical.T.Stage) + factor(Clinical.N.Stage), data = lungOpenMale, family="binomial")
summary(openDensityM)

openDensityF <- glm(nintyDay~ sm_area + age + factor(clinical.T.Stage) + factor(Clinical.N.Stage), data = lungOpenFemale, family="binomial")
summary(openDensityF)



openDensitySS <- glm(nintyDay~ densityMedianSS + age + factor(clinical.T.Stage) + factor(Clinical.N.Stage), data = lungOpenNew, family="binomial")
summary(openDensitySS)

openAreaSS <- glm(nintyDay~ densityAreaSS + age + factor(clinical.T.Stage) + factor(Clinical.N.Stage), data = lungOpenNew, family="binomial")
summary(openAreaSS)

lungOpenNew_tt <- lungOpenNew%>%
  #filter(Clinical.M.Stage == 0)
  filter(clinical.T.Stage < 5) %>%
  filter(Clinical.N.Stage < 4) 

  
sarcSurv <- coxph(Surv(time = Survival.time, event = deadstatus.event)~densityAreaSS + age + factor(gender) + factor(clinical.T.Stage) + factor(Clinical.N.Stage), data = lungOpenNew_tt)
summary(sarcSurv)

sarcSurv2 <- coxph(Surv(time = Survival.time, event = deadstatus.event)~densityAreaSS, data = lungOpenNew_tt)
sarcSurv2 <- coxph(Surv(time = Survival.time, event = deadstatus.event)~densityMedianSS, data = lungOpenNew_tt)
summary(sarcSurv2)


