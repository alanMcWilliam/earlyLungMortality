
### 28th Sept updated to git repository 
### 18th Oct added samity check image removal - spinal level and segmentation

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(summarytools)
library(tidyr)


######################################################################


lung20 <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/20_fraction_pt_vr2_chemo_vr2.csv')
lung20halfMask <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/alan_half_mask_stats.csv')
lung20fullMask <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/alan_stats.csv')

locationSanity <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/20fracLocationSanity.csv')
segSanity <- read.csv('C:/Users/alan_/Desktop/Lung_sarcopenia/20fracSanity.csv')

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

### sanity check - remove images

lung20 <- lung20 %>%
  filter(!ID %in% locationSanity$ID) %>%
  filter(!ID %in% segSanity$ID)


#### sumarise data and plot some graphs

summary(lung20$SM.Area)
summary(lung20$SM.Density)

ggplot(data=lung20, aes(halfDensity)) + 
  geom_histogram(breaks=seq(-20,40, by = 2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = "", x = "Muscle density" ) +
  theme(panel.background = element_blank())

ggplot(data=lung20, aes(fullDensity)) + 
  geom_histogram(breaks=seq(-20,40, by = 2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = "", x = "Muscle density" ) +
  theme(panel.background = element_blank())

ggplot(data = lung20) +
  geom_histogram(aes(x = halfDensity), 
                 alpha=0.3, fill ="red",binwidth=2,position="dodge") +
  geom_histogram(aes(x = fullDensity), 
                 alpha=0.3, fill ="green",binwidth=2,position="dodge") +
  labs(title = "", x = "Muscle density") +
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
  mutate(T.stageClean = case_when(T.stage == 'T1' ~ '1',
                                  T.stage == 'T1a' ~ '1',
                                  T.stage == 'T1b' ~ '1',
                                  T.stage == 'T2' ~ '2',
                                  T.stage == 'T2a' ~ '2', 
                                  T.stage == 'T2b' ~ '2',
                                  T.stage == 'T3' ~ '3',
                                  T.stage == 'T4' ~ '4')) %>%
  mutate(N.stageClean = case_when(N.stage == 'N0' ~ '0',
                                  N.stage == 'N1' ~ '1', 
                                  N.stage == 'N2' ~ '2', 
                                  N.stage == 'N2c' ~ '2',
                                  N.stage == 'N3' ~ '3')) %>%
  mutate_at(vars(nintyDay, hundredEighty), ~replace_na(.,FALSE)) %>%
  select(Age, tumour.size, fullArea, halfDensity, nintyDay, hundredEighty, performance.status, gender, T.stageClean, N.stageClean)

view(dfSummary(lung20clean))

#lung20Pat <- lung20 %>%
#  select(ID)

lung20clean <- lung20clean[complete.cases(lung20clean),]

### convert to cm3 assuming sampled at ~1mm
lung20clean$fullAreaCC <- lung20clean$fullArea*0.01
### remove single PS 4
#lung20clean <- lung20clean[as.numeric(as.character(lung20clean$performance.status)) < 4, ]

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


sarcDensity <- glm(lung20clean$nintyDay~lung20clean$halfDensity)
summary(sarcDensity)
#sarcDensity <- glm(nintyDay~halfDensity + Age + log(tumour.size) + gender, data = lung20clean, family=binomial(link='logit'))
#summary(sarcDensity)
sarcDensity <- glm(nintyDay~halfDensity + Age + T.stageClean + N.stageClean + gender, data = lung20clean, family=binomial(link='logit'))
summary(sarcDensity)


test <- glm(lung20clean$hundredEighty~lung20clean$halfDensity)
summary(test)

summary(lung20clean$fullAreaCC)
ggplot(data=lung20clean, aes(x = fullAreaCC, color = nintyDay)) + 
  geom_density() +
  labs(title = "", x = "Muscle area" ) +
  theme(panel.background = element_blank())

wilcox.test(lung20clean$fullAreaCC~lung20clean$nintyDay)

sarcNinty <- glm(nintyDay~fullAreaCC + Age + log(tumour.size) + gender, data = lung20clean, family=binomial(link='logit'))
summary(sarcNinty)
sarcNinty <- glm(nintyDay~fullAreaCC + Age + T.stageClean + N.stageClean + gender, data = lung20clean, family=binomial(link='logit'))
summary(sarcNinty)

test2 <- glm(hundredEighty~fullAreaCC + Age + log(tumour.size) + gender, data = lung20clean, family=binomial(link='logit'))
summary(test2)

m1 <- AIC(sarcDensity)
m2 <- AIC(sarcNinty)
anova(sarcDensity, sarcNinty, test = "Chisq")
anova(nintyDayPS, sarcNinty, test = "Chisq")


summary(lung20clean$Age)
lung20clean$ageBin <- cut(lung20clean$Age, c(40, 50, 60, 70, 80, 90))

ggplot(lung20clean) +
  geom_boxplot(aes(x = ageBin, y = fullAreaCC, fill = nintyDay))

nintyDayPS <- glm(nintyDay~performance.status + Age + T.stageClean + N.stageClean  + gender, data = lung20clean, family=binomial(link='logit'))
summary(nintyDayPS)

p <-ggplot(lung20clean) +
  geom_boxplot(aes(x = performance.status, y = fullAreaCC, fill = nintyDay))
show(p)

ggplot(lung20clean) +
  geom_boxplot(aes(x = gender, y = fullAreaCC, fill = nintyDay)) +
  labs(title = "", y = "Muscle area" ) +
  theme(panel.background = element_blank())

ggplot(lung20clean) +
  geom_boxplot(aes(x = gender, y = halfDensity, fill = nintyDay)) +
  labs(title = "", y = "Muscle density" ) +
  theme(panel.background = element_blank())

tapply(lung20clean$nintyDay, lung20clean$gender, summary)
tapply(lung20clean$nintyDay, lung20clean$gender, summary)

t <- lung20clean %>%
  filter(gender == "Male")
wilcox.test(t$fullAreaCC~t$nintyDay)
wilcox.test(t$halfDensity~t$nintyDay)


ggplot(t) +
  geom_boxplot(aes(x = nintyDay, y = fullAreaCC, fill = nintyDay))


t <- lung20clean %>%
  filter(gender == "Female")
wilcox.test(t$fullAreaCC~t$nintyDay)
wilcox.test(t$halfDensity~t$nintyDay)

+
  stat_compare_means(aes(group = gender))


p <- ggboxplot(lung20clean, x = "gender", y = "fullAreaCC",
               color = "nintyDay", palette = "jco",
               add = "jitter")
p + stat_compare_means(aes(group = nintyDay))


###################################################
#### try and do some cross-validation

library(boot)

cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
model <- glm(nintyDay~fullAreaCC + Age + T.stageClean + N.stageClean + gender, data = lung20clean, family=binomial(link='logit'))

t <- 1-cv.glm(lung20clean, model,K=5,cost=cost)$delta[1]
t

lungClean2 <- lung20clean %>%
  select(Age, fullAreaCC, nintyDay, gender,T.stageClean, N.stageClean, fullAreaCC)
lungClean2b <- lung20clean %>%
  select(Age, fullAreaCC, nintyDay, gender,T.stageClean, N.stageClean, halfDensity)

# Define training control
set.seed(123) 
train.control <- trainControl(method = "cv", number = 5, savePredictions = T)
# Train the model
model <- train(as.factor(nintyDay)~., data = lungClean2b,
               trControl = train.control, "rf", preProc=c("center", "scale"))
# Summarize the results
print(model)
#plot(model)

### nothing here works...
library(pROC)
selectedIndices <- model$pred$mtry == 2
plot.roc(model$pred$obs[selectedIndices],
         model$pred$M[selectedIndices])

library(MLeval)
res <- evalm(model)
res$roc


library(plotROC)
g <- ggplot(model$pred[selectedIndices, ], aes(m=M, d=factor(obs, levels = c("R", "M")))) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc()
plot(g)
###################################################
library(randomForest)


lungClean3 <- lung20clean %>%
  select(Age, fullAreaCC, nintyDay, gender,T.stageClean, N.stageClean, fullAreaCC, performance.status, halfDensity)

sarcRF <- randomForest(as.factor(nintyDay) ~ fullAreaCC + Age + T.stageClean + N.stageClean + gender + performance.status, ntree=500, data=lungClean3, importance=TRUE,
                        proximity=TRUE)

sarcRF
importance(sarcRF, type = 1)



#######################
#I've put sanity check images for the 55/20 cohort on server 2 (/data/Donal/lung_55_20_SANITY)
#I'd already gone through and checked the level labelling so that's what location_sanity
#+ location_sanity.csv is. I had tried two ways of getting the coordinates so DSNT just means pick the second option. (All of this should be accounted for when calculating stats) (edited) 
#For those where the level labelling landed on a different vert. I estimated where L3 might be by adding/removing an offset - so some segmentations are bound to be at the wrong level
#Let me know if/when you want me to do this for Fabio's data
