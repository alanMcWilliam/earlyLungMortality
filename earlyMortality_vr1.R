
### 28th Sept updated to git repository 


library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(summarytools)


######################################################################


lung20 <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/20_fraction_pt_vr2_chemo_vr2.csv')
lung20halfMask <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/alan_half_mask_stats.csv')
lung20fullMask <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia//alan_stats.csv')



######################################################################

lung20halfMask$ID <- str_remove(lung20halfMask$ID, '.npy')
lung20halfMask <- lung20halfMask %>%
  rename(halfArea = "SM.Area",
         halfDensity = "SM.Density")

lung20fullMask$ID <- str_remove(lung20fullMask$ID, '.npy')
lung20fullMask <- lung20fullMask %>%
  rename(fullArea = "SM.Area",
         fullDensity = "SM.Density")


lung20 <- merge(lung20, lung20fullMask, by = 'ID')
lung20 <- merge(lung20, lung20halfMask, by = 'ID')
View(lung20)

#### sumarise data and plot some graphs

summary(lung20$SM.Area)
summary(lung20$SM.Density)

ggplot(data=lung20, aes(halfDensity)) + 
  geom_histogram(breaks=seq(-20,40, by = 2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "Muscle density" ) +
  theme(panel.background = element_blank())

ggplot(data=lung20, aes(fullDensity)) + 
  geom_histogram(breaks=seq(-20,40, by = 2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "Muscle density" ) +
  theme(panel.background = element_blank())

ggplot(data = lung20) +
  geom_histogram(aes(x = halfDensity), 
                 alpha=0.3, fill ="red",binwidth=2,position="dodge") +
  geom_histogram(aes(x = fullDensity), 
                 alpha=0.3, fill ="green",binwidth=2,position="dodge") +
  labs(title = i, x = "Muscle density") +
  theme(panel.background = element_blank())

plot(lung20$fullDensity, lung20$halfDensity)


###########################################################################

lung20$dead_at2 <- as.numeric(as.character(lung20$dead_at))
summary(lung20$dead_at2)

lung20$nintyDay <- as.numeric(as.character(lung20$dead_at)) <= 3
lung20$hundredEighty <- as.numeric(as.character(lung20$dead_at)) <= 6

summary(lung20$nintyDay)
summary(lung20$hundredEighty)

###clean up dataset
lung20clean <- lung20 %>%
  select(Age, tumour.size, fullArea, halfDensity, nintyDay, hundredEighty, performance.status, gender, T.stage, N.stage) %>%
  mutate(T.stageClean = case_when(T.stage == 'T1' ~ 'T1',
                                  T.stage == 'T1a' ~ 'T1',
                                  T.stage == 'T1b' ~ 'T1',
                                  T.stage == 'T2' ~ 'T2',
                                  T.stage == 'T2a' ~ 'T2', 
                                  T.stage == 'T2b' ~ 'T2',
                                  T.stage == 'T3' ~ 'T3',
                                  T.stage == 'T4' ~ 'T4')) %>%
  mutate(N.stageClean = case_when(N.stage == 'N0' ~ 'N0',
                                  N.stage == 'N1' ~ 'N1', 
                                  N.stage == 'N2' ~ 'N2', 
                                  N.stage == 'N2c' ~ 'N2',
                                  N.stage == 'N3' ~ 'N3')) %>%
  select(Age, tumour.size, fullArea, halfDensity, nintyDay, hundredEighty, performance.status, gender, T.stageClean, N.stageClean)

lung20clean$comp <- complete.cases(lung20clean)
lung20clean <- lung20clean[lung20clean$comp == TRUE,]
lung20clean <- lung20clean[-ncol(lung20clean)]

### convert to cm3 assuming sampled at ~1mm
lung20clean$fullAreaCC <- lung20clean$fullArea*0.01
### remove single PS 4
lung20clean <- lung20clean[as.numeric(as.character(lung20clean$performance.status)) < 4, ]

view(dfSummary(lung20clean))

ggplot(data=lung20clean, aes(x = halfDensity, color = nintyDay)) + 
  geom_histogram(breaks=seq(-20,40, by = 2)) +
  labs(title = "", x = "Muscle density" ) +
  theme(panel.background = element_blank())

ggplot(data=lung20clean, aes(x = halfDensity, color = hundredEighty)) + 
  geom_histogram(breaks=seq(-20,40, by = 2)) +
  labs(title = "", x = "Muscle density" ) +
  theme(panel.background = element_blank())

wilcox.test(lung20clean$halfDensity~lung20clean$nintyDay)
wilcox.test(lung20clean$halfDensity~lung20clean$hundredEighty)


test <- glm(lung20clean$nintyDay~lung20clean$halfDensity)
summary(test)
test <- glm(lung20clean$hundredEighty~lung20clean$halfDensity)
summary(test)

summary(lung20clean$fullAreaCC)
ggplot(data=lung20clean, aes(x = fullAreaCC, color = nintyDay)) + 
  geom_histogram(breaks=seq(0,100, by = 5)) +
  labs(title = "", x = "Muscle area" ) +
  theme(panel.background = element_blank())

wilcox.test(lung20clean$fullAreaCC~lung20clean$nintyDay)
sarcNinty <- glm(nintyDay~fullAreaCC + Age + log(tumour.size) + gender, data = lung20clean, family=binomial(link='logit'))
summary(sarcNinty)
sarcNinty <- glm(nintyDay~fullAreaCC + Age + T.stageClean + N.stageClean + gender, data = lung20clean, family=binomial(link='logit'))
summary(sarcNinty)

test2 <- glm(hundredEighty~fullAreaCC + Age + tumour.size, data = lung20clean, family=binomial(link='logit'))
summary(test2)


summary(lung20clean$Age)
lung20clean$ageBin <- cut(lung20clean$Age, c(45, 50, 55, 60, 65, 70, 75, 80, 85, 90))

ggplot(lung20clean) +
  geom_boxplot(aes(x = ageBin, y = fullAreaCC, fill = nintyDay))

nintyDayPS <- glm(nintyDay~performance.status + Age + T.stageClean + N.stageClean  + gender, data = lung20clean, family=binomial(link='logit'))
summary(nintyDayPS)

ggplot(lung20clean) +
  geom_boxplot(aes(x = performance.status, y = fullAreaCC, fill = nintyDay))





#######################
#I've put sanity check images for the 55/20 cohort on server 2 (/data/Donal/lung_55_20_SANITY)
#I'd already gone through and checked the level labelling so that's what location_sanity
#+ location_sanity.csv is. I had tried two ways of getting the coordinates so DSNT just means pick the second option. (All of this should be accounted for when calculating stats) (edited) 
#For those where the level labelling landed on a different vert. I estimated where L3 might be by adding/removing an offset - so some segmentations are bound to be at the wrong level
#Let me know if/when you want me to do this for Fabio's data
