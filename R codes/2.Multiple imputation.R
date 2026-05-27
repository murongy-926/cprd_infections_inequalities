# Multiple imputation 

library(bit64)
library(dplyr)
library(data.table)
library(splines)
library(ggplot2)
library(beepr)

install.packages("mice")
install.packages("future")
install.packages("rlang")
install.packages("furrr")
library(furrr)
library(mice)
library(future)
library(furrr)



################################################################################
# 1. Data preparation 
################################################################################

# Open dataset 

master <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_1619_with_missing_covar.rds")

summary(is.na(master))

#------------------------------------------------------------------------------#
# (1)    Generate auxiliary variable - comorbidity 
#------------------------------------------------------------------------------#

# Select comorbidity <= end of 2018

comorbid <- fread("D:\\CPRDData\\Analysis_Nam\\Final datasets 19Jan\\Datasets to run codes\\0-visit-comorbid data.txt")

comed <- fread("D:\\CPRDData\\Analysis_Nam\\Final datasets 19Jan\\Datasets to run codes\\0-visit-comed data.txt")

# Select comorbidities 

table(comorbid$group)

comorbid_select <- comorbid %>%
  filter(group %in% c("diab", "heart", "liv", "neu", "obes", "renal", "res"))

comorbid_select <- comorbid_select %>%
  mutate(obsdate = as.Date(obsdate)) %>%
  filter(obsdate < as.Date("2018-12-31"))

# Select patients with imd 

table(comed$group)

comed_select <- comed %>%
  filter(group == "imd")

comed_select <- comed_select %>%
  mutate(obsdate = as.Date(obsdate)) %>%
  filter(obsdate < as.Date("2018-12-31"))

# Combine datasets

comorbid_all <- rbind(comorbid_select, comed_select)

# Generate one row per patient (patients with comorbidities that should be excluded from analysis)

comorbid_pat <- comorbid_all %>%
  group_by(patid) %>%
  slice(1)

# Save dataset

saveRDS(comorbid_pat, "D:\\CPRDData\\Analysis_Murong\\Data\\Multiple imputation\\comorbid_pat.rds") #N=2627341

# Add comorbidity to master dataset

summary(is.na(comorbid_pat$group)) #no missing

master <- master %>%
  left_join(comorbid_pat %>% select(patid, group), by = "patid")

summary(is.na(master$group)) 

master <- master %>%
  mutate(comorbid = case_when(
    !is.na(group) ~ "Yes",
    TRUE ~ "No"
  )) %>%
  rename(comorbid_group = group) 

table(master$comorbid)

master$comorbid <- factor(master$comorbid, levels = c("No", "Yes"))


#------------------------------------------------------------------------------#
# (2)    Clean smoking status "M" category
#------------------------------------------------------------------------------#

master <- master %>%
  mutate(current_smkstatus_update = case_when(
    current_smkstatus_update == "M" ~ NA_character_,
    TRUE ~ current_smkstatus_update
  ))

table(master$current_smkstatus_update)

master$current_smkstatus_update <- factor(
  master$current_smkstatus_update,
  levels = c("non-smoker", "ex-smoker", "current smoker")
)

saveRDS(master, "D:/CPRDData/Analysis_Murong/Data/master_1619_with_missing_covar.rds")

#------------------------------------------------------------------------------#
# (2)    Separate child and adult 
#------------------------------------------------------------------------------#

# Check variable class
sapply(master, class)

# Update variable class

master$vacc_flu <- factor(master$vacc_flu, 
                          levels = c(0, 1), 
                          labels = c("No", "Yes"))


master$vacc_pneu <- factor(master$vacc_pneu, 
                           levels = c(0, 1), 
                           labels = c("No", "Yes"))

# Generate child & adult datasets 

child <- master %>% filter(age <18)
adult <- master %>% filter(age >=18)

bs_age <- bs(child$age, df = 4)
colnames(bs_age) <- paste0("bs_age_", seq_len(ncol(bs_age)))
attr_age <- attributes(bs_age)
child <- cbind(child, bs_age)


bs_age <- bs(adult$age, df = 4)
colnames(bs_age) <- paste0("bs_age_", seq_len(ncol(bs_age)))
attr_age <- attributes(bs_age)
adult <- cbind(adult, bs_age)


#------------------------------------------------------------------------------#
# (3)    Select variables for MI
#------------------------------------------------------------------------------#

child_mi <- child %>%
  select(patid, gender, bs_age_1, bs_age_2, bs_age_3, bs_age_4, ethnic_2, imd, vacc_pneu, vacc_flu, comorbid)

summary(is.na(child_mi))

adult_mi <- adult %>%
  select(patid, gender, bs_age_1, bs_age_2, bs_age_3, bs_age_4, ethnic_2, imd, vacc_pneu, vacc_flu, bmi, current_smkstatus_update, comorbid)

summary(is.na(adult_mi))

saveRDS(child_mi, "D:\\CPRDData\\Analysis_Murong\\Data\\Multiple imputation\\child_for_mi.rds") 
saveRDS(adult_mi, "D:\\CPRDData\\Analysis_Murong\\Data\\Multiple imputation\\adult_for_mi.rds") 


################################################################################
# 2. MI for children 
################################################################################

# Check class of variables
sapply(child_mi, class)

# Keep only the variables used in the imputation model
imp_dat <- child_mi[, c("imd", "ethnic_2", "gender",
                  "bs_age_1", "bs_age_2", "bs_age_3", "bs_age_4",
                  "vacc_pneu", "vacc_flu", "comorbid")]

# Set up mice defaults
ini  <- mice(imp_dat, maxit = 0)

meth <- ini$method
pred <- ini$predictorMatrix

# Model for each variable  
meth[] <- ""
meth["imd"] <- "pmm"
meth["ethnic_2"] <- "logreg"

# Predict imd and ethnic_2 using all the other variables
pred[,] <- 0
pred["imd", ] <- c(
  imd = 0,
  ethnic_2 = 1,
  gender = 1,
  bs_age_1 = 1,
  bs_age_2 = 1,
  bs_age_3 = 1,
  bs_age_4 = 1,
  vacc_pneu = 1,
  vacc_flu = 1,
  comorbid = 1
)

pred["ethnic_2", ] <- c(
  imd = 1,
  ethnic_2 = 0,
  gender = 1,
  bs_age_1 = 1,
  bs_age_2 = 1,
  bs_age_3 = 1,
  bs_age_4 = 1,
  vacc_pneu = 1,
  vacc_flu = 1,
  comorbid = 1
)


Start_time <- Sys.time()

imp <- futuremice(
  imp_dat,
  m = 20,
  maxit = 10,
  method = meth,
  predictorMatrix = pred,
  parallelseed = 123,
  n.core = 20
)

End_time <- Sys.time()

# Quick diagnostics
plot(imp, c("imd", "ethnic_2"))
densityplot(imp, ~ imd)
stripplot(imp, imd ~ .imp, pch = 20, cex = 0.6)
bwplot(imp, imd ~ .imp)


completed_all <- complete(imp, "broad")

imd_cols <- grep("^imd\\.", names(completed_all), value = TRUE)
eth_cols <- grep("^ethnic_2\\.", names(completed_all), value = TRUE)

completed_subset <- completed_all[, c(imd_cols, eth_cols)]

# Obtain imputed values across completed datasets

set.seed(123)

completed_subset$imd_imp <- apply(
  completed_subset[, imd_cols],
  1,
  function(x) sample(x, 1)
)


completed_subset <- completed_subset %>%
  mutate(ethnic_2_imp = ifelse(
           rowMeans(across(starts_with("ethnic_2."), ~ as.numeric(. == "Ethnic minority"))) > 0.5,
           1, 0
         )
         ) %>%
  mutate(ethnic_2_imp = factor(ethnic_2_imp, levels = c(0, 1), labels = c("White", "Ethnic minority"))) %>%
  mutate(imd_imp = as.integer(imd_imp))


# Generate clean imputed dataset (only keep variables for subsequent analysis)
raw_dataset_select <- child %>%
  select(patid, imd, ethnic_2, gender, 
         age, bs_age_1, bs_age_2, bs_age_3, bs_age_4,
         vacc_pneu, vacc_flu, comorbid)


imputed_dataset <- cbind(raw_dataset_select, completed_subset %>% select(imd_imp, ethnic_2_imp))

sapply(imputed_dataset, class)

# check IMD distribution

ggplot(imputed_dataset) +
  geom_histogram(aes(x = imd, fill = "Raw"),
                 alpha = 0.5, position = "identity", bins = 30, na.rm = TRUE) +
  geom_histogram(aes(x = imd_imp, fill = "Imputed"),
                 alpha = 0.5, position = "identity", bins = 30, na.rm = TRUE) +
  labs(title = "Distribution of IMD",
       x = "IMD",
       y = "Count",
       fill = "Legend") +
  theme_minimal()

saveRDS(imp, "D:/CPRDData/Analysis_Murong/Data/Multiple imputation/child_imputed_mids.rds")
saveRDS(imputed_dataset, "D:/CPRDData/Analysis_Murong/Data/Multiple imputation/child_imputed_clean.rds")




################################################################################
# 3. MI for adults 
################################################################################

sapply(adult_mi, class)

# Keep only the variables used in the imputation model
imp_dat <- adult_mi[, c("imd", "ethnic_2", "gender",
                        "bs_age_1", "bs_age_2", "bs_age_3", "bs_age_4",
                        "vacc_pneu", "vacc_flu", "bmi", "current_smkstatus_update", "comorbid")]

summary(is.na(imp_dat)) #missing: imd, ethnic_2, bmi, current_smkstatus_update

# Set up mice defaults
ini  <- mice(imp_dat, maxit = 0)

meth <- ini$method
pred <- ini$predictorMatrix

# Define model for variables
meth[] <- ""
meth["imd"] <- "pmm"
meth["ethnic_2"] <- "logreg"
meth["bmi"] <- "pmm"
meth["current_smkstatus_update"] <- "polyreg"    


# Predict imd ethnic_2 smoking bmi using all the other variables
pred[,] <- 0
pred["imd", ] <- c(
  imd = 0,
  ethnic_2 = 1,
  gender = 1,
  bs_age_1 = 1,
  bs_age_2 = 1,
  bs_age_3 = 1,
  bs_age_4 = 1,
  vacc_pneu = 1,
  vacc_flu = 1,
  current_smkstatus_update = 1,
  bmi = 1,
  comorbid = 1
)

pred["ethnic_2", ] <- c(
  imd = 1,
  ethnic_2 = 0,
  gender = 1,
  bs_age_1 = 1,
  bs_age_2 = 1,
  bs_age_3 = 1,
  bs_age_4 = 1,
  vacc_pneu = 1,
  vacc_flu = 1,
  current_smkstatus_update = 1,
  bmi = 1,
  comorbid = 1
)

pred["bmi", ] <- c(
  imd = 1,
  ethnic_2 = 1,
  gender = 1,
  bs_age_1 = 1,
  bs_age_2 = 1,
  bs_age_3 = 1,
  bs_age_4 = 1,
  vacc_pneu = 1,
  vacc_flu = 1,
  current_smkstatus_update = 1,
  bmi = 0,
  comorbid = 1
)

pred["current_smkstatus_update", ] <- c(
  imd = 1,
  ethnic_2 = 1,
  gender = 1,
  bs_age_1 = 1,
  bs_age_2 = 1,
  bs_age_3 = 1,
  bs_age_4 = 1,
  vacc_pneu = 1,
  vacc_flu = 1,
  current_smkstatus_update = 0,
  bmi = 1,
  comorbid = 1
)


# Run multiple imputation

install.packages("mice")
install.packages("future")
install.packages("rlang")
install.packages("furrr")
library(furrr)
library(mice)
library(future)



Start_time <- Sys.time()

imp <- futuremice(
  imp_dat,
  m = 20,
  maxit = 10,
  method = meth,
  predictorMatrix = pred,
  parallelseed = 123,
  n.core = 20
)

End_time <- Sys.time()

End_time-Start_time #2.2 hours

beep("ping")


# Quick diagnostics
plot(imp, c("imd", "ethnic_2"))
#stripplot(imp, imd ~ .imp, pch = 20, cex = 0.6)
densityplot(imp, ~ imd)
bwplot(imp, imd ~ .imp)

# Distribution looks very similar between datasets, check whether they are identical
all.equal(completed_all$imd.1, completed_all$imd.2) #"Mean relative difference: 0.643176", not identical 



# Completed datasets

completed_all <- complete(imp, "broad")

imd_cols <- grep("^imd\\.", names(completed_all), value = TRUE)
eth_cols <- grep("^ethnic_2\\.", names(completed_all), value = TRUE)
bmi_cols <- grep("^bmi\\.", names(completed_all), value = TRUE)
smk_cols <- grep("^current_smkstatus_update\\.", names(completed_all), value = TRUE)

completed_subset <- completed_all[, c(imd_cols, eth_cols, bmi_cols, smk_cols)]

summary(completed_subset)

# IMD 

set.seed(123)

completed_subset$imd_imp <- apply(
  completed_subset[, imd_cols],
  1,
  function(x) sample(x, 1)
)


completed_subset <- completed_subset %>%
  mutate(ethnic_2_imp = ifelse(
           rowMeans(across(starts_with("ethnic_2."), ~ as.numeric(. == "Ethnic minority"))) > 0.5,
           1, 0
         ),
         bmi_imp = round(rowMeans(across(starts_with("bmi."))), 2)
  ) %>%
  mutate(ethnic_2_imp = factor(ethnic_2_imp, levels = c(0, 1), labels = c("White", "Ethnic minority")))%>%
  mutate(imd_imp = as.integer(imd_imp))



# For smoking 

smk_cols <- grep("^current_smkstatus_update\\.", names(completed_subset), value = TRUE)

completed_subset$current_smkstatus_update_imp <- apply(
  completed_subset[, smk_cols],
  1,
  function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
)

completed_subset$current_smkstatus_update_imp <- factor(
  completed_subset$current_smkstatus_update_imp,
  levels = c("non-smoker", "ex-smoker", "current smoker")
)

raw_dataset_select <- adult %>%
  select(patid, imd, ethnic_2, gender, 
         age, bs_age_1, bs_age_2, bs_age_3, bs_age_4,
         vacc_pneu, vacc_flu, bmi, current_smkstatus_update, comorbid)

imputed_dataset <- cbind(raw_dataset_select, completed_subset %>% select(imd_imp, ethnic_2_imp, bmi_imp, current_smkstatus_update_imp))

sapply(imputed_dataset, class)

# check IMD distribution

ggplot(imputed_dataset) +
  geom_histogram(aes(x = imd, fill = "Raw"),
                 alpha = 0.5, position = "identity", bins = 30, na.rm = TRUE) +
  geom_histogram(aes(x = imd_imp, fill = "Imputed"),
                 alpha = 0.5, position = "identity", bins = 30, na.rm = TRUE) +
  labs(title = "Distribution of IMD",
       x = "IMD",
       y = "Count",
       fill = "Legend") +
  theme_minimal()


saveRDS(imp, "D:/CPRDData/Analysis_Murong/Data/Multiple imputation/adult_imputed_mids.rds")
saveRDS(imputed_dataset, "D:/CPRDData/Analysis_Murong/Data/Multiple imputation/adult_imputed_clean.rds")


