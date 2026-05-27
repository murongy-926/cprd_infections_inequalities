# Estimate number of missing non-consulting patients and estimated total patients CPRD

library(bit64)
library(dplyr)
library(beepr)
library(ggplot2)
library(parallel)
library(doParallel)

options(scipen = 999)

################################################################################
# Expand national data for single age 95+ based on CPRD 95+ distribution
################################################################################

# 1. Open datasets ###

# Number of England GP registered patients by imd*single age*sex*ethnicity in July 2019

df_single <- readRDS("D:/CPRDData/Analysis_Murong/Data/National data/gp_imd_age_single_sex_eth_fitted.rds") 

sum(df_single$n_patients) 

# Master dataset with consulting patients only 

master <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_1619_with_missing_covar.rds")


# 2. Master dataset: number of patients by imd*single age*sex*ethnicity in July 2019

# Number of patients present on July 2019

master$midyear_present1 <- ifelse((master$startdate<=as.Date("2019-07-01") & master$enddate>as.Date("2019-07-01")),1,0)
summary(master$midyear_present1)  
table(master$midyear_present1) # 

summary(is.na(master$midyear_present1))

# Keep mid-year (July 1 2019) population

cprd <- master %>%
  filter(midyear_present1 == 1) %>%
  dplyr::select(patid, pracid, imd, gender, age, ethnic_2) %>%
  rename(sex = gender) %>%
  mutate()

# 3. Compare cprd vs national data

str(cprd)
str(df_single)

# Make variable characters consistent

df_single <- df_single %>%
  rename(age = age_num) %>%
  mutate(imd = as.integer(imd),
         age = as.numeric(age),
         sex = as.factor(sex), 
         ethnic_2 = as.factor(ethnic_2))

cprd <- cprd %>%
  mutate(imd = as.integer(imd),
         age = as.numeric(age),
         sex = as.factor(sex), 
         ethnic_2 = as.factor(ethnic_2))


# Number of patients by imd*single age*sex*ethnicity

cprd_grid <- cprd %>%
  group_by(imd, age, sex,ethnic_2) %>%
  count()   


# Check elements in each category 

summary(cprd_grid) 
summary(df_single) 

# 4. Expand df_single to age >95 years

cprd_95plus <- cprd_grid %>%
  filter(age >= 95) %>%
  group_by(age, imd, sex, ethnic_2) %>%
  summarise(n_cprd = sum(n), .groups = "drop") 
                                            

nat_95plus <- df_single %>%
  filter(age == 95) %>%
  select(-age) %>%
  rename(n_nat = n_patients)

sum(nat_95plus$n_nat) 


# Calculate proportion by age*sex in CPRD data (instead of by age*sex*imd*ethnicity which is ideally but may have zeros)

cprd_95plus_age_sex <- cprd_95plus %>%
  group_by(age, sex) %>%
  summarise(n_cprd_sex_age = sum(n_cprd), .groups = "drop") %>%
  group_by(sex) %>%
  mutate(prop_age_within_sex = n_cprd_sex_age / sum(n_cprd_sex_age)) %>%
  ungroup()

sum(cprd_95plus$n_cprd) 
sum(cprd_95plus_age_sex$n_cprd_sex_age) 

check_props <- cprd_95plus_age_sex %>% group_by(sex) %>% summarise(sum_prop = sum(prop_age_within_sex))



# Generate national grid data adding cprd 95+ age

imd_levels <- unique(df_single$imd)
sex_levels <- unique(df_single$sex)
ethnic_levels <- unique(df_single$ethnic_2)
cprd_old_ages <- sort(unique(cprd$age[cprd$age >= 95]))

nat_95plus_age_sex <- expand.grid(
  age = cprd_old_ages,
  sex = sex_levels,
  imd = imd_levels,
  ethnic_2 = ethnic_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) 

nat_95plus_age_sex <- nat_95plus_age_sex %>%
  left_join(nat_95plus, by = c("imd", "sex", "ethnic_2"))

sum(nat_95plus_age_sex$n_nat) 

# --- join CPRD age-by-sex proportions to the national grid ------------------

# join by age and sex (these props are the same for all IMD/ethnic combos)
nat_95plus_age_sex <- nat_95plus_age_sex %>%
  left_join(cprd_95plus_age_sex %>% select(age, sex, prop_age_within_sex),
            by = c("age", "sex"))

# --- allocate national 95+ totals into single ages ---------------------------

nat_95plus_age_sex <- nat_95plus_age_sex %>%
  mutate(n_nat_age = n_nat * prop_age_within_sex)

summary(is.na(nat_95plus_age_sex$n_nat_age))


nat_95plus_age_sex <- nat_95plus_age_sex %>%
  mutate(n_nat_age = case_when(
    !is.na(n_nat_age) ~ n_nat_age,
    is.na(n_nat_age) ~ 0
  ))

sum(nat_95plus_age_sex$n_nat_age) 


# Expand national data

df_single_1 <- df_single %>%
  filter(age < 95)

df_single_2 <- nat_95plus_age_sex %>%
  rename(n_patients = n_nat_age) %>%
  select(imd, age, sex, ethnic_2, n_patients)

df_nat <- bind_rows(df_single_1, df_single_2)

sum(df_nat$n_patients) 

saveRDS(df_nat, "D:/CPRDData/Analysis_Murong/Data/National data/gp_imd_age_single_expand_sex_eth_fitted_update.rds")


################################################################################
# Estimate missing consulting patients in CPRD 
################################################################################

# 1. Calculate number of estimated cprd registered patients

practice <- read.csv("D:/CPRDData/Analysis_Murong/Data/National data/gp-reg-pat-prac-lsoa-all.csv") 

# Number of practices 

n_gp <- practice %>%
  group_by(PRACTICE_CODE) %>%
  count() 

nat_prac <- 6904 


# Average number of registered patients by imd*single age*sex*ethnicity per practice

df_nat <- df_nat %>%
  mutate(avg_n_patients = n_patients/ nat_prac)

summary(df_nat$avg_n_patients)

# Estimated number of registered patients by imd*single age*sex*ethnicity in CPRD 

# Number of practices in cprd 

cprd_prac <- master %>%
  group_by(pracid) %>%
  count()

cprd_prac <- 1686

df_nat <- df_nat %>%
  mutate(n_cprd_est_pat = round(avg_n_patients * cprd_prac, 0))

sum(df_nat$n_cprd_est_pat) #14617881


# 2. Calculate number of consulting patients in CPRD master dataset

df_cprd <- master %>%
  filter(midyear_present1 == 1) %>%
  select(patid, age, gender, ethnic_2, imd) %>%
  group_by(imd, age, gender, ethnic_2) %>%
  count() %>%
  rename(n_cprd_consult_pat = n) %>%
  ungroup() 



# 3. Calculate number of non-consulting patients in cprd 
df_merge <- df_nat%>%
  rename(gender = sex) %>%
  left_join(df_cprd, by = c("imd", "age", "gender", "ethnic_2"))
  
# check missing 
summary(is.na(df_merge)) 

# replace missing with 0 
df_merge <- df_merge %>%
  mutate(n_cprd_consult_pat = ifelse(is.na(n_cprd_consult_pat), 0, n_cprd_consult_pat))

# Calculate missing cprd patients
df_merge <- df_merge %>%
  # Number of non-consulting patients in each imd*single age*sex*ethnicity
  mutate(n_cprd_miss_pat = n_cprd_est_pat - n_cprd_consult_pat)

# check
sum(df_merge$n_cprd_miss_pat) 
sum(df_merge$n_cprd_miss_pat)/sum(df_merge$n_cprd_est_pat)*100 

summary(df_merge$n_cprd_miss_pat)

check <- df_merge %>%
  filter(n_cprd_miss_pat <0)

sum(check$n_cprd_miss_pat) 

417/sum(df_merge$n_cprd_est_pat)*100 

# rename and save dataset
df_merge <- df_merge %>%
  rename(n_nat_pat = n_patients)

saveRDS(df_merge, "D:/CPRDData/Analysis_Murong/Data/National data/nat_cprd_non_consult_patients_update.rds")




################################################################################
# Add non-consulting patients to master dataset 
################################################################################

# 1. Generate aggregate dataset 

# (1) Assign weights (for gamlss)

df_agg <- df_merge %>%
  select(imd, age, gender, ethnic_2, n_cprd_miss_pat) %>%
  rename(weights = n_cprd_miss_pat) %>%
  # replace negative weights with 0
  mutate(weights = ifelse(weights <0, 0, weights)) 

summary(df_agg$weights) 
sum(df_agg$weights) 

# Exclude subgroup with 0 weights - no patients are missing from CPRD in this subgroup

df_agg <- df_agg %>%
  filter(weights >0) 


# (2) Generate fake patid

max_real_id <- max(as.numeric(master$patid), na.rm = TRUE)
max_real_id 

df_agg <- df_agg %>%
  arrange(imd, age, gender, ethnic_2) %>%
  mutate(patid = as.integer64(99000000000000+row_number())) %>% 
  select(patid, everything())


# (3) Add exposure, outcome, time at risk, age group into the dataset

# IMD quintile

df_agg$imd_quintile <- ifelse(df_agg$imd <= 4, 1,
                           ifelse(df_agg$imd <= 8, 2,
                                  ifelse(df_agg$imd <= 12, 3,
                                         ifelse(df_agg$imd <= 16, 4,
                                                5))))

df_agg <- df_agg %>%
  mutate(imd_quintile=recode(imd_quintile,
                             `1`="Quintile 1(least deprived)",
                             `2`="Quintile 2",
                             `3`="Quintile 3",
                             `4`="Quintile 4",
                             `5`="Quintile 5"))


# outcome & time at risk 

df_agg <- df_agg %>%
  mutate(
    # Outcome
    n_new_case = 0,
    n_new_abx = 0,
    abyes = 0,    
    
    # Time at risk
    time_at_risk = 365,
    time_at_risk_abx = 365,
    log_timeatrisk = log(time_at_risk),
    log_timeatrisk_abx = log(time_at_risk_abx)

  )

# age group - for incidence calculation

assign_age_group <- function(age) {
  if (age >=0 & age <=1) {
    return("<1")
  } else if (age >= 1 & age <= 4) {
    return("1-4")
  } else if (age >= 5 & age <= 9) {
    return("5-9")
  } else if (age >= 10 & age <= 14) {
    return("10-14")
  } else if (age >= 15 & age <= 19) {
    return("15-19")
  } else if (age >= 20 & age <= 24) {
    return("20-24")
  } else if (age >= 25 & age <= 29) {
    return("25-29")
  } else if (age >= 30 & age <= 34) {
    return("30-34")
  } else if (age >= 35 & age <= 39) {
    return("35-39")
  } else if (age >= 40 & age <= 44) {
    return("40-44")
  } else if (age >= 45 & age <= 49) {
    return("45-49")
  } else if (age >= 50 & age <= 54) {
    return("50-54")
  } else if (age >= 55 & age <= 59) {
    return("55-59")
  } else if (age >= 60 & age <= 64) {
    return("60-64")
  } else if (age >= 65 & age <= 69) {
    return("65-69")
  } else if (age >= 70 & age <= 74) {
    return("70-74")
  } else if (age >= 75 & age <= 79) {
    return("75-79")
  } else if (age >= 80 & age <= 84) {
    return("80-84")
  } else if (age >= 85 & age <= 89) {
    return("85-89")
  } else if (age >= 90 & age <= 94) {
    return("90-94")
  } else if (age >= 95) {
    return("95+")
  } else {
    return(NA)
  }
}


assign_age_group_uk <- function(age) {
  if (age >= 0 & age <= 4) {
    return("0-4")
  } else if (age >= 5 & age <= 9) {
    return("5-9")
  } else if (age >= 10 & age <= 14) {
    return("10-14")
  } else if (age >= 15 & age <= 17) {
    return("15-17")
  } else if (age >= 18 & age <= 24) {
    return("18-24")
  } else if (age >= 25 & age <= 29) {
    return("25-29")
  } else if (age >= 30 & age <= 34) {
    return("30-34")
  } else if (age >= 35 & age <= 39) {
    return("35-39")
  } else if (age >= 40 & age <= 44) {
    return("40-44")
  } else if (age >= 45 & age <= 49) {
    return("45-49")
  } else if (age >= 50 & age <= 54) {
    return("50-54")
  } else if (age >= 55 & age <= 59) {
    return("55-59")
  } else if (age >= 60 & age <= 64) {
    return("60-64")
  } else if (age >= 65 & age <= 69) {
    return("65-69")
  } else if (age >= 70 & age <= 74) {
    return("70-74")
  } else if (age >= 75 & age <= 79) {
    return("75-79")
  } else if (age >= 80 & age <= 84) {
    return("80-84")
  } else if (age >= 85 & age <= 89) {
    return("85-89")
  } else if (age >= 90) {
    return("90+")
  } else {
    return(NA)
  }
}


df_agg <- df_agg %>%
  mutate(age_group=sapply(age, assign_age_group)) %>%
  mutate(age_group_uk=sapply(age, assign_age_group_uk))

# save dataset

saveRDS(df_agg, "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_strata_update.rds") 
                                                                             
#------------------------------------------------------------------------------#

# (4) Add contributors - Option 1: use the level from master dataset 

# Method 1: Hot-deck (row sampling) -use this for GAMLSS

# ----------- Duplicate dataset to one row per individual  ----------------

df_agg_duplicate <- df_agg %>%
  select(-c("patid")) %>%
  uncount(weights) # obs=7768564

df_agg_duplicate <- df_agg_duplicate %>%
  mutate(patid = as.integer64(99000000000000+row_number())) %>%
  select(patid, everything())

names(df_agg_duplicate)

# generate child/adults data

df_agg_duplicate_child <- df_agg_duplicate %>%
  filter(age < 18)

df_agg_duplicate_adult<- df_agg_duplicate %>%
  filter(age >= 18)

# save dataset
saveRDS(df_agg_duplicate, "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds") 




################################################################################
# Generate BMI, smoking, vaccination for non-consulting patients               #
################################################################################

# Use England data


library(bit64)
library(dplyr)
library(gtsummary)
library(parallel)
library(beepr)

#------------------------------------------------------------------------------#
#                                Smoking                                       #
#------------------------------------------------------------------------------#

smk_pre_final <- readRDS("D:/CPRDData/Analysis_Murong/Data/National data/Contributor/smk_pre_final.rds")
cohort <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds") #N=7768564

cohort <- cohort %>% 
  filter(age >=18) 

table(cohort$age_group_uk)
table(smk_pre_final$age)

# Generate IMD decile

cohort$imd_decile <- ceiling(cohort$imd / 2)

library(dplyr)
set.seed(123)

cohort_imp <- cohort %>%
  mutate(across(c(imd_decile, age_group_uk, gender, ethnic_2), as.character)) %>%
  left_join(smk_pre_final %>% 
              rename(imd_decile = imd,
                     age_group_uk = age,
                     gender = sex,
                     ethnic_2 = eth) %>%
              mutate(across(c(imd_decile, age_group_uk, gender, ethnic_2), as.character)) %>%
              dplyr::select(-c(cur_smokers, ex_smokers)),
            by = c("imd_decile","age_group_uk","gender","ethnic_2")) 


cohort_imp <- cohort_imp %>%
  mutate(u = runif(n()),     #randomisation here 
         smoke_imp = case_when(
           is.na(p_cur_smoker) ~ NA_character_,
           u <= p_cur_smoker ~ "current smoker",
           u <= p_cur_smoker + p_ex_smoker ~ "ex-smoker",
           TRUE ~ "non-smoker"
         )) 




# Check 

check <- cohort_imp %>%
  select(imd_quintile, age, gender, ethnic_2, smoke_imp) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",         # Percent only for categorical
      all_continuous() ~ "{mean} ({sd})"  # Mean & SD for continuous
    ))

# Rename smoking variable to be consistent with master dataset

cohort_imp <- cohort_imp %>%
  rename(current_smkstatus_update = smoke_imp)

# Check missing
summary(is.na(cohort_imp$current_smkstatus_update)) #no missing

# Save dataset

saveRDS(cohort_imp, "D:/CPRDData/Analysis_Murong/Data/Covariates/Smoking/non_consult_smk_england_update.rds")

# Save clean dataset

smk_clean <- cohort_imp %>%
  dplyr::select(patid, current_smkstatus_update)

saveRDS(smk_clean, "D:/CPRDData/Analysis_Murong/Data/Covariates/Smoking/non_consult_smk_england_clean_update.rds")

# Merge it to master dataset




#------------------------------------------------------------------------------#
#                        Pneumococcal vaccine                                  #
#------------------------------------------------------------------------------#

pneu_pre_final <- readRDS("D:/CPRDData/Analysis_Murong/Data/National data/Contributor/pneu_pre_final.rds")
cohort <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds")

table(cohort$age_group_uk)
table(pneu_pre_final$age)

# Generate age group used in pneu prevalence dataset

cohort <- cohort %>%
  mutate(age_group_pneu = case_when(
    
    age <= 69 ~ as.character(age),
    age >= 70 & age <=74 ~ "70-74",
    age >=75 ~ "75+",
    TRUE ~ NA_character_
    
  ))


# Generate IMD decile

cohort$imd_decile <- ceiling(cohort$imd / 2)

cohort_imp <- cohort %>%
  mutate(across(c(imd_decile, age_group_pneu), as.character)) %>%
  left_join(pneu_pre_final %>% 
              rename(age_group_pneu = age) %>%
              dplyr::select(imd_decile, age_group_pneu, prevalence) %>%
              mutate(across(c(imd_decile, age_group_pneu), as.character)),
            by = c("imd_decile","age_group_pneu")) 

summary(is.na(cohort_imp$prevalence)) 

cohort_imp <- cohort_imp %>%
  mutate(prevalence = ifelse(is.na(prevalence), 0, prevalence))

set.seed(123)

cohort_imp <- cohort_imp %>%
  mutate(
    vacc_pneu = rbinom(n(), size = 1, prob = prevalence),
    vacc_pneu = if_else(vacc_pneu == 1, "Yes", "No")
  )

table(cohort_imp$vacc_pneu)


# Check 

check <- cohort_imp %>%
  filter(age >=18) %>%  
  select(imd_quintile, vacc_pneu) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",         # Percent only for categorical
      all_continuous() ~ "{mean} ({sd})"  # Mean & SD for continuous
    ))


# Check missing
summary(is.na(cohort_imp$vacc_pneu)) #no missing

# Save dataset

saveRDS(cohort_imp, "D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/non_consult_pneu_england_update.rds")

# Save clean dataset

pneu_clean <- cohort_imp %>%
  dplyr::select(patid, age, vacc_pneu)

saveRDS(pneu_clean, "D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/non_consult_pneu_england_clean_update.rds")

# Merge it to master dataset





#------------------------------------------------------------------------------#
#                                       BMI                                    #
#------------------------------------------------------------------------------#

bmi_final <- readRDS("D:/CPRDData/Analysis_Murong/Data/National data/Contributor/bmi_final.rds")
cohort <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds") 

table(cohort$age_group)
table(bmi_final$age)

# Generate age group used in pneu prevalence dataset

cohort <- cohort %>%
  filter(age >= 18) %>%    
  mutate(age_group_bmi = case_when(
    
    age >=18 & age <=19 ~ "18-19",
    age >=20 & age <=24 ~ "20-24",
    age >=25 & age <=29 ~ "25-29",
    age >=30 & age <=34 ~ "30-34",
    age >=35 & age <=39 ~ "35-39",
    age >=40 & age <=44 ~ "40-44",
    age >=45 & age <=49 ~ "45-49",
    age >=50 & age <=54 ~ "50-54",
    age >=55 & age <=59 ~ "55-59",
    age >=60 & age <=64 ~ "60-64",
    age >=65 & age <=69 ~ "65-69",
    age >=70 & age <=74 ~ "70-74",
    age >=75 & age <=79 ~ "75-79",
    age >=80 & age <=84 ~ "80-84",
    age >=85 & age <=89 ~ "85-89",
    age >=90  ~ "90+",
    TRUE ~ NA_character_
    
  ))


# Generate IMD decile

cohort$imd_decile <- ceiling(cohort$imd / 2)

# Check number of strata

check <- cohort %>%
  group_by(imd_decile, age_group_bmi, gender, ethnic_2) %>%
  count() #obs=638, while bmi has 314 strata


# ---------- Prepare BMI summaries at different levels

bmi <- readRDS("D:/CPRDData/Analysis_Murong/Data/National data/Contributor/bmi_from_stata.rds")

# Label imd to be consistent with cprd data  
bmi <- bmi %>%
  mutate(imd = case_when(
    
    imd == 1 ~ "Quintile 1(least deprived)",
    imd == 2 ~ "Quintile 2",
    imd == 3 ~ "Quintile 3",
    imd == 4 ~ "Quintile 4",
    imd == 5 ~ "Quintile 5"
    
  ))


bmi_full <- bmi %>%
  group_by(imd, age, sex, eth) %>%
  summarise(
    bmi_mean = mean(bmi, na.rm = TRUE),
    bmi_sd   = sd(bmi, na.rm = TRUE),
    n        = n(),
    .groups = "drop"
  )

bmi_imd_age_sex <- bmi %>%
  group_by(imd, age, sex) %>%
  summarise(
    bmi_mean_i_a_s = mean(bmi, na.rm = TRUE),
    bmi_sd_i_a_s   = sd(bmi, na.rm = TRUE),
    .groups = "drop"
  )

bmi_imd_age_eth <- bmi %>%
  group_by(imd, age, eth) %>%
  summarise(
    bmi_mean_i_a_e = mean(bmi, na.rm = TRUE),
    bmi_sd_i_a_e   = sd(bmi, na.rm = TRUE),
    .groups = "drop"
  )

bmi_age_sex <- bmi %>%
  group_by(age, sex) %>%
  summarise(
    bmi_mean_a_s = mean(bmi, na.rm = TRUE),
    bmi_sd_a_s   = sd(bmi, na.rm = TRUE),
    .groups = "drop"
  )

bmi_age <- bmi %>%
  group_by(age) %>%
  summarise(
    bmi_mean_a = mean(bmi, na.rm = TRUE),
    bmi_sd_a   = sd(bmi, na.rm = TRUE),
    .groups = "drop"
  )


# ------  Join all summaries to cohort_imp

cohort_imp <- cohort %>%
  dplyr::select(patid, imd_quintile, imd_decile, age_group_bmi, gender, ethnic_2) %>%
  mutate(across(c(imd_quintile, age_group_bmi, gender, ethnic_2), as.character)) %>%
  #rename(bmi_cprd = bmi) %>%
  left_join(bmi_full %>%
              mutate(across(c(imd, age, sex, eth), as.character)),
            by = c("imd_quintile" = "imd",
                   "age_group_bmi" = "age",
                   "gender" = "sex",
                   "ethnic_2" = "eth")) %>%
  left_join(bmi_imd_age_sex %>%
              mutate(across(c(imd, age, sex), as.character)),
            by = c("imd_quintile" = "imd",
                   "age_group_bmi" = "age",
                   "gender" = "sex")) %>%
  left_join(bmi_imd_age_eth %>%
              mutate(across(c(imd, age, eth), as.character)),
            by = c("imd_quintile" = "imd",
                   "age_group_bmi" = "age",
                   "ethnic_2" = "eth")) %>%
  left_join(bmi_age_sex %>%
              mutate(across(c(age, sex), as.character)),
            by = c("age_group_bmi" = "age",
                   "gender" = "sex")) %>%
  left_join(bmi_age %>%
              mutate(across(c(age), as.character)),
            by = c("age_group_bmi" = "age"))


summary(is.na(cohort_imp))


# -----  Coalesce to pick the most granular available value

cohort_imp <- cohort_imp %>%
  mutate(
    bmi_mean_final = coalesce(
      bmi_mean,
      bmi_mean_i_a_s,
      bmi_mean_i_a_e,
      bmi_mean_a_s,
      bmi_mean_a
    ),
    bmi_sd_final = coalesce(
      bmi_sd,
      bmi_sd_i_a_s,
      bmi_sd_i_a_e,
      bmi_sd_a_s,
      bmi_sd_a
    ),
    bmi_source = case_when(
      !is.na(bmi_mean)        ~ "imd_age_sex_eth",
      !is.na(bmi_mean_i_a_s)  ~ "imd_age_sex",
      !is.na(bmi_mean_i_a_e)  ~ "imd_age_eth",
      !is.na(bmi_mean_a_s)    ~ "age_sex",
      !is.na(bmi_mean_a)      ~ "age",
      TRUE                    ~ "missing"
    )
  )

table(is.na(cohort_imp$bmi_mean_final)) #no missing
summary(cohort_imp$bmi_mean_final)

table(cohort_imp$bmi_source)

#------------------------------------------------------------------------------#
# Estimate BMI value 

library(truncnorm)

set.seed(123)
cohort_imp <- cohort_imp %>%
  mutate(
    bmi = rtruncnorm(
      n(),
      a = 15,
      b = 60,
      mean = bmi_mean_final,
      sd   = bmi_sd_final
    )
  )
summary(cohort_imp$bmi) 

summary(cohort_imp$bmi_sd_final) 

cohort_imp <- cohort_imp %>%
  mutate(
    bmi_sd_adj = pmin(bmi_sd_final, 8),  
    bmi_cap = rtruncnorm(n(),
                         a = 15,
                         b = 60,
                         mean = bmi_mean_final,
                         sd = bmi_sd_adj)
  )


summary(cohort_imp$bmi_cap) 


# Check 

check <- cohort_imp %>%
  select(imd_quintile, bmi, bmi_cap, bmi_mean_final) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",         # Percent only for categorical
      all_continuous() ~ "{mean} ({sd})"  # Mean & SD for continuous
    ))

check # almost no difference between bmi and bmi_cap among imd quintile, use bmi 

check <- bmi %>%
  group_by(imd) %>%
  summarise(bmi_mean = mean(bmi), .groups = "drop")


# Check missing
summary(is.na(cohort_imp$bmi)) #no missing

# Save dataset

saveRDS(cohort_imp, "D:/CPRDData/Analysis_Murong/Data/Covariates/BMI/non_consult_bmi_england_update.rds")

# Save clean dataset

bmi_clean <- cohort_imp %>%
  dplyr::select(patid, bmi)

saveRDS(bmi_clean, "D:/CPRDData/Analysis_Murong/Data/Covariates/BMI/non_consult_bmi_england_clean_update.rds")

# Merge it to master dataset


#------------------------------------------------------------------------------#
#                        Flu vaccine                                           #
#------------------------------------------------------------------------------#

flu_pre_final <- readRDS("D:/CPRDData/Analysis_Murong/Data/National data/Contributor/flu_pre_final.rds")
cohort <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds")


# Generate age group used in flu prevalence dataset (only for 90+)

cohort <- cohort %>%
  mutate(age_group_flu = case_when(
    
    age < 90 ~ as.character(age),
    age >=90 ~ "90+",
    TRUE ~ NA_character_
    
  ))


flu_pre_final <- flu_pre_final %>%
  mutate(
    age_group_flu = case_when(
      as.numeric(as.character(age)) < 90  ~ as.character(age),
      as.numeric(as.character(age)) == 90 ~ "90+",
      TRUE ~ NA_character_
    )
  )


# Generate IMD decile

cohort$imd_decile <- ceiling(cohort$imd / 2)

cohort_imp <- cohort %>%
  mutate(across(c(imd_decile, age_group_flu), as.character)) %>%
  left_join(flu_pre_final %>% 
              dplyr::select(imd_decile, age_group_flu, prevalence) %>%
              mutate(across(c(imd_decile, age_group_flu), as.character)),
            by = c("imd_decile","age_group_flu")) 

summary(is.na(cohort_imp$prevalence)) # no missing


set.seed(123)

cohort_imp <- cohort_imp %>%
  mutate(
    vacc_flu = rbinom(n(), size = 1, prob = prevalence),
    vacc_flu = if_else(vacc_flu == 1, "Yes", "No")
  )

table(cohort_imp$vacc_flu)


# Check 

check <- cohort_imp %>%
  filter(age < 18) %>%  
  select(imd_quintile,vacc_flu_cprd, vacc_flu) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",         # Percent only for categorical
      all_continuous() ~ "{mean} ({sd})"  # Mean & SD for continuous
    ))


# Check missing
summary(is.na(cohort_imp$vacc_flu)) #no missing

# Save dataset

saveRDS(cohort_imp, "D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/non_consult_flu_england_update.rds")

# Save clean dataset

flu_clean <- cohort_imp %>%
  dplyr::select(patid, age, vacc_flu)

saveRDS(flu_clean, "D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/non_consult_flu_england_clean_update.rds")

# Merge it to master dataset


#-------------------------------------------------------------------------------
# Add contributors to non-consulting master dataset 
#-------------------------------------------------------------------------------

df_agg_duplicate <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds") #N=7768564

names(df_agg_duplicate)


# Covariates dataset
smk <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/Smoking/non_consult_smk_england_clean_update.rds")
flu <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/non_consult_flu_england_clean_update.rds")
bmi <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/BMI/non_consult_bmi_england_clean_update.rds")
pneu <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/non_consult_pneu_england_clean_update.rds")

df_agg_duplicate <- df_agg_duplicate %>%
  # Smoking
  left_join(smk, by = "patid") %>%
  # BMI 
  left_join(bmi, by = "patid") %>%
  # Pneu
  left_join(pneu %>% 
              dplyr::select(-age), by = "patid") %>%
  # Flu
  left_join(flu %>% 
              dplyr::select(-age), by = "patid") %>%
  # Set to factor variable
  mutate(current_smkstatus_update = factor(
    current_smkstatus_update,
    levels = c("non-smoker", "ex-smoker", "current smoker") #non-smoker as reference
  )) %>%
  mutate(vacc_flu = factor(
    vacc_flu,
    levels = c("No", "Yes") #No as reference
  )) %>%
  mutate(vacc_pneu = factor(
    vacc_pneu,
    levels = c("No", "Yes") #No as reference
  )) 


# Save dataset  
saveRDS(df_agg_duplicate, "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds") #N=7768564


################################################################################
# Add non-consult-patient dataset to master & all disease master dataset
################################################################################

# 1.Master dataset
master_imputed <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds") 


df_agg_duplicate <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds") 

names(df_agg_duplicate)
names(master_imputed)

common_vars <- intersect(names(master_imputed), names(df_agg_duplicate))
common_vars

master_added <- bind_rows(master_imputed %>% dplyr::select(all_of(common_vars)) %>% mutate(source = "consult"), 
                          df_agg_duplicate %>% dplyr::select(all_of(common_vars)) %>% mutate(source = "non_consult"))


table(master_added$source) 

saveRDS(master_added,"D:/CPRDData/Analysis_Murong/Data/master_consult_nonconsult_pat_with_england_contributors_update.rds")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# 2.  CHILDREN disease master

# ---------- Adding missing patients to all disease master dataset -------------

# (1) Function for both-sex diseases

add_missing_patient <- function(infection){
  
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    
  consult_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/"
  add_data_path <- "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds"  #THIS INCLUDES BOTH CHILDREN AND ADULTS 
  save_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/"
  
  # Upload original master data
  consult_data <- readRDS(paste0(consult_data_path, infection, "_master_child.rds"))
  assign("consult_data", consult_data, envir = .GlobalEnv)
  
  on.exit({
    rm(consult_data, envir = .GlobalEnv)
    gc()
  }, add = TRUE)
  
  
  # Upload added data
  add_data <- readRDS(paste0(add_data_path))
  assign("add_data", add_data, envir = .GlobalEnv)
  
  on.exit({
    rm(add_data, envir = .GlobalEnv)
    gc()
  }, add = TRUE)
  

  logs$steps <- c(logs$steps, "Uploaded data")
  
  # bind datasets
  
  common_vars <- intersect(names(consult_data), names(add_data))
  
  
  full_data <- bind_rows(consult_data %>% dplyr::select(all_of(common_vars)) %>% mutate(source = "consult"), 
                         add_data %>% filter(age <18) %>% dplyr::select(all_of(common_vars)) %>% mutate(source = "non_consult")
                         )
  
  logs$steps <- c(logs$steps, "Binded data")
  
  # save merged dataset
  
  saveRDS(full_data, paste0(save_data_path, infection, "_master_child.rds"))
  
  list(success = TRUE, logs = logs)
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))
    
    
  })
    
}
    

infections <- c("asthma", "urti", "lrti", "cough", "pneu", "cell", "impe", "exter", "media", "gast",  
                "sore", "sinu", "uurti", "bron", "ulrti", "copd") #copd and pros only used for calculating incidence rate, exclude utio, pyel and pros as they are one-sex, their syntax see below


num_cores <- 20

cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(dplyr)
})


clusterExport(cl, c("add_missing_patient", "infections"))

#Run parallel  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- parLapply(cl, infections, add_missing_patient) 

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")



# (2) Function for one-sex diseases

# Female
add_missing_patient_female <- function(infection){
  
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    
    consult_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/"
    add_data_path <- "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds"  #THIS INCLUDES BOTH CHILDREN AND ADULTS 
    save_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/"
    
    # Upload original master data
    consult_data <- readRDS(paste0(consult_data_path, infection, "_master_child.rds"))
    assign("consult_data", consult_data, envir = .GlobalEnv)
    
    on.exit({
      rm(consult_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    # Upload added data
    add_data <- readRDS(paste0(add_data_path))
    assign("add_data", add_data, envir = .GlobalEnv)
    
    on.exit({
      rm(add_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # bind datasets
    
    common_vars <- intersect(names(consult_data), names(add_data))
    
    
    full_data <- bind_rows(consult_data %>% 
                             dplyr::select(all_of(common_vars)) %>% mutate(source = "consult"), 
                           
                           add_data %>% 
                             filter(age <18) %>%
                             filter(gender == "Female") %>% #UPDATE HERE 
                             dplyr::select(all_of(common_vars)) %>% mutate(source = "non_consult")
    )
  
    
    logs$steps <- c(logs$steps, "Binded data")
    
    # save merged dataset
    
    saveRDS(full_data, paste0(save_data_path, infection, "_master_child.rds"))
    
    list(success = TRUE, logs = logs)
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- rlang::trace_back()
    return(list(success = FALSE, logs = logs))
    
    
  })
  
}


# Male 

add_missing_patient_male <- function(infection){
  
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    
    consult_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/"
    add_data_path <- "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds"  #THIS INCLUDES BOTH CHILDREN AND ADULTS 
    save_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/"
    
    # Upload original master data
    consult_data <- readRDS(paste0(consult_data_path, infection, "_master_child.rds"))
    assign("consult_data", consult_data, envir = .GlobalEnv)
    
    on.exit({
      rm(consult_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    # Upload added data
    add_data <- readRDS(paste0(add_data_path))
    assign("add_data", add_data, envir = .GlobalEnv)
    
    on.exit({
      rm(add_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # bind datasets
    
    common_vars <- intersect(names(consult_data), names(add_data))
    
    
    full_data <- bind_rows(consult_data %>% 
                             dplyr::select(all_of(common_vars)) %>% mutate(source = "consult"), 
                           
                           add_data %>% 
                             filter(age < 18) %>%
                             filter(gender == "Male") %>% #UPDATE HERE 
                             dplyr::select(all_of(common_vars))%>% mutate(source = "non_consult") 
    )
    
    
    logs$steps <- c(logs$steps, "Binded data")
    
    # save merged dataset
    
    saveRDS(full_data, paste0(save_data_path, infection, "_master_child.rds"))
    
    list(success = TRUE, logs = logs)
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- rlang::trace_back()
    return(list(success = FALSE, logs = logs))
    
    
  })
  
}


# Run function 

infections <- c("utio", "pyel")

results <- lapply(infections, add_missing_patient_female) 

beep("ping")

infections <- c("pros")

results <- lapply(infections, add_missing_patient_male) 

beep("ping")


# ----- Generate spline terms & log time at risk datasets for GAMLSS -----------

generate_var_child <- function(infection) {
  
  # Define input and output paths to avoid overwriting original data
  input_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/", infection, "_master_child.rds") 
  
  # Create an output directory if it doesn't exist
  output_dir <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Data/Child/"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  
  infection_data <- readRDS(input_path)
  
  # Generate B-spline for imd
  bs_imd <- bs(infection_data$imd, df = 2, degree = 2) 
  colnames(bs_imd) <- paste0("bs_imd_", seq_len(ncol(bs_imd)))
  attr_imd <- attributes(bs_imd)
  
  # Generate B-spline for age
  bs_age <- bs(infection_data$age, df = 4)
  colnames(bs_age) <- paste0("bs_age_", seq_len(ncol(bs_age)))
  attr_age <- attributes(bs_age)
  
  # Combine with original data
  infection_data <- cbind(infection_data, bs_imd, bs_age)
  
  # Generate log time at risk 
  
  infection_data <- infection_data %>%
    mutate(log_timeatrisk = log(time_at_risk),
           log_timeatrisk_abx = log(time_at_risk_abx))
  
  # Save to the new, safe location
  output_path <- paste0(output_dir, infection, "_master_child.rds") 
  output_path_imd <- paste0(output_dir, infection, "_child_attr_imd.rds") 
  output_path_age <- paste0(output_dir, infection, "_child_attr_age.rds") 
  
  saveRDS(infection_data, output_path)
  saveRDS(attr_imd, output_path_imd)
  saveRDS(attr_age, output_path_age)
  
  return(TRUE) # Explicitly return TRUE for success
}


# --- Parallel Setup ---

num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Ensure the cluster is always stopped
on.exit(stopCluster(cl))


# --- Run in Parallel ---


infections <- c("asthma", "urti", "lrti", "cough", "pneu", "cell", "impe", "exter", "media", "gast", "utio", "pyel", 
                "sore", "sinu", "uurti", "bron", "ulrti")
results <- foreach(
  infection = infections,
  .packages = c("splines", "dplyr"),      
  .errorhandling = "pass"     
) %dopar% {
  generate_var_child(infection)   
}

print("Processing summary:")
print(setNames(results, infections))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# 3. ADULTS disease master  

# ------------ Adding missing patients to all disease master datasets ---------- 

# (1) Both-sex diseases

add_missing_patient <- function(infection){
  
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    
    consult_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/"
    add_data_path <- "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds"  #THIS INCLUDES BOTH CHILDREN AND ADULTS 
    save_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/"
    
    # Upload original master data
    consult_data <- readRDS(paste0(consult_data_path, infection, "_master_adult.rds"))
    assign("consult_data", consult_data, envir = .GlobalEnv)
    
    on.exit({
      rm(consult_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    # Upload added data
    add_data <- readRDS(paste0(add_data_path))
    assign("add_data", add_data, envir = .GlobalEnv)
    
    on.exit({
      rm(add_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # bind datasets
    
    common_vars <- intersect(names(consult_data), names(add_data))
    
    
    full_data <- bind_rows(consult_data %>% dplyr::select(all_of(common_vars)), 
                           add_data %>% filter(age >=18) %>% dplyr::select(all_of(common_vars))
    )
    
    logs$steps <- c(logs$steps, "Binded data")
    
    # save merged dataset
    
    saveRDS(full_data, paste0(save_data_path, infection, "_master_adult.rds"))
    
    list(success = TRUE, logs = logs)
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))
    
    
  })
  
}


infections <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti") 

num_cores <- 20

cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(dplyr)
})


clusterExport(cl, c("add_missing_patient", "infections"))

# Run parallel  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- parLapply(cl, infections, add_missing_patient) 

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")



# (2) One-sex diseases

add_missing_patient_female <- function(infection){
  
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    
    consult_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/"
    add_data_path <- "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds"  #THIS INCLUDES BOTH CHILDREN AND ADULTS 
    save_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/"
    
    # Upload original master data
    consult_data <- readRDS(paste0(consult_data_path, infection, "_master_adult.rds"))
    assign("consult_data", consult_data, envir = .GlobalEnv)
    
    on.exit({
      rm(consult_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    # Upload added data
    add_data <- readRDS(paste0(add_data_path))
    assign("add_data", add_data, envir = .GlobalEnv)
    
    on.exit({
      rm(add_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # bind datasets
    
    common_vars <- intersect(names(consult_data), names(add_data))
    
    
    full_data <- bind_rows(consult_data %>% 
                             dplyr::select(all_of(common_vars)), 
                           
                           add_data %>% 
                             filter(age >= 18) %>%
                             filter(gender == "Female") %>% #UPDATE HERE 
                             dplyr::select(all_of(common_vars)) 
    )
    
    
    logs$steps <- c(logs$steps, "Binded data")
    
    # save merged dataset
    
    saveRDS(full_data, paste0(save_data_path, infection, "_master_adult.rds"))
    
    list(success = TRUE, logs = logs)
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- rlang::trace_back()
    return(list(success = FALSE, logs = logs))
    
    
  })
  
}



add_missing_patient_male <- function(infection){
  
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    
    consult_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/"
    add_data_path <- "D:/CPRDData/Analysis_Murong/Data/master_non_consult_pat_individual_update.rds"  #THIS INCLUDES BOTH CHILDREN AND ADULTS 
    save_data_path <- "D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/"
    
    # Upload original master data
    consult_data <- readRDS(paste0(consult_data_path, infection, "_master_adult.rds"))
    assign("consult_data", consult_data, envir = .GlobalEnv)
    
    on.exit({
      rm(consult_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    # Upload added data
    add_data <- readRDS(paste0(add_data_path))
    assign("add_data", add_data, envir = .GlobalEnv)
    
    on.exit({
      rm(add_data, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # bind datasets
    
    common_vars <- intersect(names(consult_data), names(add_data))
    
    
    full_data <- bind_rows(consult_data %>% 
                             dplyr::select(all_of(common_vars)), 
                           
                           add_data %>% 
                             filter(age >= 18) %>%
                             filter(gender == "Male") %>% #UPDATE HERE 
                             dplyr::select(all_of(common_vars)) 
    )
    
    
    logs$steps <- c(logs$steps, "Binded data")
    
    # save merged dataset
    
    saveRDS(full_data, paste0(save_data_path, infection, "_master_adult.rds"))
    
    list(success = TRUE, logs = logs)
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- rlang::trace_back()
    return(list(success = FALSE, logs = logs))
    
    
  })
  
}



infections <- c("utio", "pyel")

results <- lapply(infections, add_missing_patient_female) 

beep("ping")

infections <- c("pros")

results <- lapply(infections, add_missing_patient_male) 

beep("ping")


# ----- Generate spline terms & log time at risk datasets for GAMLSS -----------

generate_var_adult <- function(infection) {
  
  # Define input and output paths to avoid overwriting original data
  input_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed_added/", infection, "_master_adult.rds") 
  
  # Create an output directory if it doesn't exist
  output_dir <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Data/Adult/"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  
  infection_data <- readRDS(input_path)
  
  # Generate B-spline for imd
  bs_imd <- bs(infection_data$imd, df = 2, degree = 2) 
  colnames(bs_imd) <- paste0("bs_imd_", seq_len(ncol(bs_imd)))
  attr_imd <- attributes(bs_imd)
  
  # Generate B-spline for age
  bs_age <- bs(infection_data$age, df = 4)
  colnames(bs_age) <- paste0("bs_age_", seq_len(ncol(bs_age)))
  attr_age <- attributes(bs_age)
  
  # Generate B-spline for bmi
  bs_bmi <- bs(infection_data$bmi, df = 5)
  colnames(bs_bmi) <- paste0("bs_bmi_", seq_len(ncol(bs_bmi)))
  attr_bmi <- attributes(bs_bmi)
  
  # Combine with original data
  infection_data <- cbind(infection_data, bs_imd, bs_age, bs_bmi)
  
  # Generate log time at risk 
  
  infection_data <- infection_data %>%
    mutate(log_timeatrisk = log(time_at_risk),
           log_timeatrisk_abx = log(time_at_risk_abx))
  
  # Save to the new, safe location
  output_path <- paste0(output_dir, infection, "_master_adult.rds") 
  output_path_imd <- paste0(output_dir, infection, "_adult_attr_imd.rds") 
  output_path_age <- paste0(output_dir, infection, "_adult_attr_age.rds") 
  output_path_bmi <- paste0(output_dir, infection, "_adult_attr_bmi.rds") 
 
  
  saveRDS(infection_data, output_path)
  saveRDS(attr_imd, output_path_imd)
  saveRDS(attr_age, output_path_age)
  saveRDS(attr_bmi, output_path_bmi)
  
  return(TRUE) # Explicitly return TRUE for success
}



# --- Parallel Setup ---

num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Ensure the cluster is always stopped
on.exit(stopCluster(cl))


# --- Run in Parallel ---

infections <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", "utio", "pyel", "pros",
                "sore", "sinu", "uurti", "bron", "ulrti")

# Note the addition of the .packages argument
results <- foreach(
  infection = infections,
  .packages = c("splines", "dplyr"),      
  .errorhandling = "pass"     
) %dopar% {
  generate_var_adult(infection)   
}

# The 'results' variable will now be a list of TRUE/FALSE values,
# making it easy to see which infections were processed successfully.
print("Processing summary:")
print(setNames(results, infections))


