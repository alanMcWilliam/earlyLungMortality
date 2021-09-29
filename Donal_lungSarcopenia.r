library(tidyverse)
library(summarytools)
library(stringr)
library(ggplot2)
library(survminer)
library(survival)
library(ranger)
library(ggfortify)
library(magrittr)
library(Hmisc)
library(tibble)
library(dplyr)
library(finalfit)
library(AMR)
##------------------- CLINICAL DATA ----------

##Clean up clinical data
cols_list <- c("patient_id", "Age.start.treatment", "primary_disease_site", "T.Stage", "N.Stage", "comorbidities", "ecog_ps", "sex", "Death.Status", "Overall.Survival.Time..Days.", "Primary.Intent")
                  #"T stage", "N stage", "M stage", "Age start treatment")
data <- read.csv("./Lung_pathology.csv")
idx <- match(cols_list, names(data))
clean_data <- tibble(data[, idx])

# Clean-up column names
clean_data %<>%
  rename(
    death_status = "Death.Status",
    survival_time = "Overall.Survival.Time..Days.",
    primary_intent = "Primary.Intent",
    age = "Age.start.treatment",
    t_stage = "T.Stage",
    n_stage = "N.Stage"
  ) %>% 
  transform(sex = as.factor(sex)) %>%
  filter(primary_disease_site=="Lung (NSCLC)") %>%
  filter(primary_intent == 'Radical') %>%
  tibble()
### --------------READ ALAN'S CLINICAL DATA
alan_data <- tibble(read.csv("./20_fraction_pt_vr2_chemo_vr2.csv"))

alan_cols <- c("ID", "Age", "T.stage", "N.stage", "co.morbity.score", "performance.status", "gender", "status", "follow_up", "chemo")
alan_idx <- match(alan_cols, names(alan_data))
clean_alan <- tibble(alan_data[, alan_idx])

clean_alan %<>%
  rename(
    death_status = "status",
    survival_time_months = "follow_up",
    age = "Age",
    patient_id = "ID",
    t_stage = "T.stage",
    n_stage = "N.stage",
    comorbidities = "co.morbity.score",
    ecog_ps = "performance.status",
    sex = "gender",
  ) %>%
  mutate(survival_time = survival_time_months * 30.417) %>%
  transform(ecog_ps = as.numeric(ecog_ps)) %>%
  tibble()
  
clean_alan %<>% 
  transform(t_stage = as.numeric(gsub(".*?([0-9]+).*", "\\1", t_stage)), 
            n_stage = as.numeric(gsub(".*?([0-9]+).*", "\\1", n_stage)),
            survival_time = as.numeric(survival_time))
# CONCATENATE ALAN'S CLINICAL DATA TO FABIO'S

clinic_data <- bind_rows(clean_data, clean_alan)

# -------------SARCOPENIA -----------
# Get sarcopenia data from analysis
sarcopenia_data <- tibble(read.csv("./stats.csv"))

sarcopenia_data_extra <- tibble(read.csv("./extra_stats.csv"))

alan_sarcopenia <- tibble(read.csv("./alan_half_mask_stats.csv"))

# Change patient id
sarcopenia_data$ID %<>% 
  str_split_fixed('.npy', 2)
# Rename cols
sarcopenia_data %<>% 
  rename(
    SMA= "SM.Area",
    SMD = "SM.Density"
  ) %>% 
  transform(SMD = as.numeric(SMD)) %>%
  tibble()

# Change patient id
sarcopenia_data_extra$ID %<>% 
  str_split_fixed('.npy', 2)
# Rename cols
sarcopenia_data_extra %<>% 
  rename(
    SMA= "SM.Area",
    SMD = "SM.Density"
  ) %>% 
  transform(SMD = as.numeric(SMD)) %>%
  tibble()

# Change patient id
alan_sarcopenia$ID %<>% 
  str_split_fixed('.npy', 2)
# Rename cols
alan_sarcopenia %<>% 
  rename(
    SMA= "SM.Area",
    SMD = "SM.Density"
  ) %>% 
  transform(SMD = as.numeric(SMD)) %>%
  tibble()


sarcopenia_data %<>% bind_rows(sarcopenia_data_extra, alan_sarcopenia)

## --------MERGE DATA
out_df <- clinic_data %>%
  merge(alan_sarcopenia, by.x='patient_id', by.y='ID') %>%
  merge(sarcopenia_data, by.x='patient_id', by.y='ID') %>%
  tibble()



## ------- CLINICAL FRAILTY SCALE

cfs_data <- tibble(read.csv("./clean_cfs.csv"))

cfs_data %<>% select(-X) %>%
  transform(cfs = as.factor(cfs)) %>% tibble()
  
# Merge frames
out_df <- clean_data %>% 
  merge(sarcopenia_data, by.x ='patient_id', by.y='ID') %>% 
  merge(cfs_data, by.x='patient_id', by.y='patient_id') %>%
  tibble()

## -----------READ IDs TO REMOVE AFTER SANITY CHECK --------
sanity<- tibble(read.csv("./seg_sanity_extra.csv"))

`%!in%` = Negate(`%in%`)
out_df %<>% filter(patient_id %!in% sanity$ID)

# Remove rows where SMD = nan
completeFun <- function(data, desiredCols){
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

out_df <- completeFun(out_df, "SMD")
out_df <- completeFun(out_df, "age")
# -------PLOTTING -----
# DATASET SUMMARY - NOT WORKING YET
summary <- dfSummary(out_df)
summarytools::view(summary)


# ----- DENSITY PLOT -----
med <- out_df %>% 
  group_by(sex) %>%
  summarise(median = median(SMD))
med


# Density split on sex
g <- ggplot(out_df, aes(x=SMD, color=sex)) + geom_density(size=2) +
  geom_vline(data=med, aes(xintercept=median, color=sex), linetype='dashed', size=1) +
  ggsave('./plots/sex_density.png')
g

# Split on age
cat_age <- out_df %>%
  mutate(age_cat = age_groups(age, split_at = c(35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95)))
males <- cat_age %>% filter(sex == "Male")
ggplot(cat_age, aes(x=age_cat, y=SMD, group=age_cat)) + 
  geom_boxplot(position=position_dodge()) +
  facet_wrap(~sex) + ggsave('./plots/smd_age.png')

# SMD vs CFS
ggplot(out_df, aes(x=cfs, y=SMD, group=cfs)) +
  geom_boxplot(position=position_dodge()) + ggsave('./plots/smd_cfs.png')

# SMD vs comorbidities
ggplot(out_df, aes(x=comorbidities, y=SMD, group=comorbidities)) +
  geom_boxplot(position=position_dodge()) + ggsave('./plots/smd_ace.png')
# SMD vs ECOG PS
ggplot(out_df, aes(x=ecog_ps, y=SMD, group=ecog_ps)) +
  geom_boxplot(position=position_dodge()) + ggsave('./plots/smd_ps.png')

ggplot(out_df, aes(x=chemo, y=SMD, group=chemo)) +
  geom_boxplot(position=position_dodge()) + ggsave('./plots/smd_tstage.png')

### -----------------KAPLAN-MEIER-----------

med_age <- out_df %>% 
  group_by(sex) %>%
  summarise(median = median(SMD))
med_age


#Categorise data on sex specific median
sex_cat_data <- out_df %>% 
  group_by(sex) %>% 
  mutate(SMD_cat = SMD > median(SMD)) %>%
  mutate(SMA_cat = SMA > median(SMA))
sex_cat_data

# Plot survival curves
sex_cat_split <- survfit(Surv(time =survival_time , event = death_status)~SMD_cat, data = sex_cat_data)
summary(sex_cat_split)
km <- ggsurvplot(sex_cat_split, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE, tables.height = 0.25,
          legend=c(0.75, 0.85)) + ggsave('./plots/KM_smd_sexMedianSplit.png')
km

# Shift females by + 1.3 such that sex-specific medians match, then split on median
females <- out_df %>%
  filter(sex =='Female') %>%
  mutate(SMD = SMD + 0.5)
males <- out_df %>% filter(sex == 'Male')
update_sex <- bind_rows(males, females)
update_sex %<>% mutate(SMD_cat = SMD > quantile(SMD, 0.5)) %>%
  mutate(quartile = ntile(SMD, 4))

cat_split <- survfit(Surv(time =survival_time , event = death_status)~SMD_cat, data = update_sex)
summary(cat_split)
km <- ggsurvplot(cat_split, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE, tables.height = 0.25,
                  legend=c(0.75, 0.85)) + ggsave('./plots/KM_smd_MedianSplit.png')
km


legend.labs= c("SMD Low", "SMD High"),
#----------------Cox proportional Hazards------------------

# UNIVARIATE
#1. SMD
out_df$SMD_cut <- out_df$SMD > 1.258
cox_smd <- coxph(Surv(time=survival_time, event=death_status)~SMD, data=out_df)
summary(cox_smd)
#2. Performance Status
cox_ps <- coxph(Surv(time=survival_time, event=death_status)~factor(ecog_ps), data=out_df)
summary(cox_ps)
#3. CFS 
cox_cfs <- coxph(Surv(time=survival_time, event=death_status)~factor(cfs), data=out_df)
summary(cox_cfs)
#4. Comorbidites
cox_com <- coxph(Surv(time=survival_time, event=death_status)~factor(comorbidities), data=out_df)
summary(cox_com)

#MULTIVARIATE
#0 Clinical variables ONLY
cox_ph_base <- coxph(Surv(time=survival_time, event=death_status) ~ sex + age + factor(t_stage) + factor(n_stage), data=out_df)
summary(cox_ph_base)
# 0.5 Add comorbidities
cox_ph_com <- coxph(Surv(time=survival_time, event=death_status) ~ factor(comorbidities) + sex + age + factor(t_stage) + factor(n_stage), data=out_df)
summary(cox_ph_com)
#1 ECOG PS
cox_ph_ps <- coxph(Surv(time=survival_time, event=death_status) ~factor(ecog_ps)  + sex + age + factor(t_stage) + factor(n_stage), data=out_df)
summary(cox_ph_ps)
#2 PS + SMD
cox_ph_smd <- coxph(Surv(time=survival_time, event=death_status) ~SMD + sex +age + factor(t_stage)+factor(n_stage), data=out_df)
summary(cox_ph_smd)
# 3 PS+ CFS
cox_ph_cfs <- coxph(Surv(time=survival_time, event=death_status) ~ sex +age + factor(t_stage)+factor(n_stage) + factor(cfs), data=out_df)
summary(cox_ph_cfs)
# 4 PS + CFS + SMD
cox_ph_all <- coxph(Surv(time=survival_time, event=death_status) ~factor(ecog_ps) + sex +age + factor(t_stage)+factor(n_stage) + SMD, data=out_df)
summary(cox_ph_all)



# Akaike information criterion (AIC)
# 0: Clinical vars, 0.5: Comorb, 1: PS , 2: PS + SMD, 3: PS + CFS, 4: PS+SMD+CFS
aic0 <- AIC(cox_ph_base)
aic05<-AIC(cox_ph_com)
aic1<-AIC(cox_ph_ps)
aic2<-AIC(cox_ph_smd)
aic3<-AIC(cox_ph_cfs)
aic4<-AIC(cox_ph_all)

###

opt_cut_mean2 <- survminer::surv_cutpoint(out_df, time = "survival_time", event = "death_status", "SMD",  minprop = 0.1, progressbar = TRUE) #smethod="logrank" set within)
summary(opt_cut_mean2)
cat_mean2 <-survminer::surv_categorize(opt_cut_mean2)
cat_split_mean2 <- survfit(Surv(time = survival_time, event = death_status)~SMD, data = cat_mean2)
ggsurvplot(cat_split_mean2, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)



### Get IDS for patients in 50-55 age group
ids55 <- cat_age %>% filter(age_cat == '50-54') %>%
  select(patient_id)

cox_age <- coxph(Surv(time=survival_time, event=death_status)~SMD + strata(age_cat), data=cat_age)
summary(cox_age)



sum(out_df$survival_time<120)

out_df$test_surv <- out_df$survival_time < 180

tapply(out_df$SMD, out_df$test_surv, summary)
tapply(out_df$SMD, out_df$test_surv, wilcox.test)
cat_age <- out_df %>%
  mutate(age_cat = age_groups(age, split_at = c(35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95)))
males <- cat_age %>% filter(sex == "Male")
ggplot(cat_age, aes(x=test_surv, y=SMD, group=test_surv)) + 
  geom_boxplot(position=position_dodge()) +
  facet_wrap(~sex) + ggsave('./plots/smd_age.png')

