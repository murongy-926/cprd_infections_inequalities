# Descriptive statistics of the study population by IMD quintile 

# Identify objects that are not named "master"
keep <- c("master", "master_child", "master_adult")

rm(list=setdiff(ls(), keep))


################################################################################
################################################################################

library(dplyr)
#install.packages("gtsummary")  
library(gtsummary)
library(gt)
library(writexl)

################################################################################
# Consulting population 
################################################################################

# Open datasets

master <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds")

# Generate age group variable 
master <- master %>%
  mutate(age_group_descriptive = case_when(
    age <= 4 ~ "0-4",
    age >= 5 & age <= 11 ~ "5-11",
    age >= 12 & age <= 17 ~ "12-17",
    age >= 18 & age <= 64 ~ "18-64",
    age >= 65 ~ "65+"
  ))


master_child <- master %>%
  filter(age < 18)

master_adult <- master %>%
  filter(age >= 18)

summary(master$age)

setwd("C:/Users/cig-murongy/Documents/Output_imputed_observed/Descriptive/")

# Children 

# age, sex, ethnicity, 
# flu, pneu, rotav, smoking, alcohol, bmi

stat_child <- master_child %>%
  select(imd_quintile, age, age_group_descriptive, gender, ethnic_2, vacc_pneu, vacc_flu) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",          # Percent only
      all_continuous()  ~ "{median} ({p25}-{p75})" # Median (IQR)
    ),
    digits = list(
      all_categorical() ~ 1,  # 1 decimal for percentages
      all_continuous()  ~ 1   # 1 decimal for mean & SD
    )
  )

    
print(stat_child)

stat_child_table <- as_tibble(stat_child)

write_xlsx(stat_child_table, "stat_child_table_observed_pop.xlsx")


# Adults 

stat_adult <- master_adult %>%
  select(imd_quintile, age, age_group_descriptive, gender, ethnic_2, vacc_pneu, vacc_flu, current_smkstatus_update, bmi) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",          # Percent only
      all_continuous()  ~ "{median} ({p25}-{p75})" # Median (IQR)
    ),
    digits = list(
      all_categorical() ~ 1,  # 1 decimal for percentages
      all_continuous()  ~ 1   # 1 decimal for mean & SD
    )
  )

print(stat_adult)

stat_adult_table <- as_tibble(stat_adult)

write_xlsx(stat_adult_table, "stat_adult_table_observed_pop.xlsx")

################################################################################
# Estimated sample 
################################################################################

# Open datasets

master_added <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_consult_nonconsult_pat_with_england_contributors_update.rds")

# Generate age group variable 
master_added <- master_added %>%
  mutate(age_group_descriptive = case_when(
    age <= 4 ~ "0-4",
    age >= 5 & age <= 11 ~ "5-11",
    age >= 12 & age <= 17 ~ "12-17",
    age >= 18 & age <= 64 ~ "18-64",
    age >= 65 ~ "65+"
  ))


master_child <- master_added %>%
  filter(age < 18)

master_adult <- master_added %>%
  filter(age >= 18)

summary(master_added$age)

setwd("C:/Users/cig-murongy/Documents/Output_imputed/Descriptive/")

# Children 

# age, sex, ethnicity, 
# flu, pneu, smoking, bmi

stat_child <- master_child %>%
  select(imd_quintile, age, age_group_descriptive, gender, ethnic_2, vacc_pneu, vacc_flu) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",          # Percent only
      all_continuous()  ~ "{median} ({p25}-{p75})" # Median (IQR)
    ),
    digits = list(
      all_categorical() ~ 1,  # 1 decimal for percentages
      all_continuous()  ~ 1   # 1 decimal for mean & SD
    )
  )

print(stat_child)

stat_child_table <- as_tibble(stat_child)

write_xlsx(stat_child_table, "stat_child_table_total_pop.xlsx")


# Adults

stat_adult <- master_adult %>%
  select(imd_quintile, age, age_group_descriptive, gender, ethnic_2, vacc_pneu, vacc_flu, current_smkstatus_update, bmi) %>%
  tbl_summary(
    by = imd_quintile,
    statistic = list(
      all_categorical() ~ "{p}",          # Percent only
      all_continuous()  ~ "{median} ({p25}-{p75})" # Median (IQR)
    ),
    digits = list(
      all_categorical() ~ 1,  # 1 decimal for percentages
      all_continuous()  ~ 1   # 1 decimal for mean & SD
    )
  )

print(stat_adult)

stat_adult_table <- as_tibble(stat_adult)

write_xlsx(stat_adult_table, "stat_adult_table_total_pop.xlsx")
