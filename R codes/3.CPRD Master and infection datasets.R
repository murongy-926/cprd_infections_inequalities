# Data preparation for analysis 
# Generate master and infections datasets for consulting patients 

setwd("D:/CPRDData/Analysis_Murong/Data")                    
dir()

install.packages("data.table")
install.packages("tidyverse")
install.packages("dplyr")

library(bit64)
library(data.table)
library(tidyverse)
library(dplyr)
library(parallel)
library(beepr)
#------------------------------------------------------------------------------#
# 1. Generate master dataset & related variables 
#------------------------------------------------------------------------------#

# Open dataset 

# Newly generated master dataset, see do file "Sample size_update.R"

master <- readRDS("D:/CPRDData/Analysis_Murong/Data/master_1619_with_missing_covar.rds") 

names(master)
summary(is.na(master))

# Check missing data: pat_present, regenddate    
# Checked, correct 


# Generate master dataset with missing covariates 

child_imputed <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/child_imputed_clean.rds")
adult_imputed <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/adult_imputed_clean.rds")

child_imputed$bmi <- NA_real_
child_imputed$current_smkstatus_update <- factor(
  NA,
  levels = levels(adult_imputed$current_smkstatus_update)
)

child_imputed$bmi_imp <- NA_real_
child_imputed$current_smkstatus_update_imp <- factor(
  NA,
  levels = levels(adult_imputed$current_smkstatus_update_imp)
)

master_imputed <- rbind(child_imputed, adult_imputed)

names(master_imputed)


master_imputed <- master_imputed %>%
  select(-c("imd", "ethnic_2", "bmi", "current_smkstatus_update")) %>%
  rename(imd = imd_imp,
         ethnic_2 = ethnic_2_imp, 
         bmi = bmi_imp,
         current_smkstatus_update = current_smkstatus_update_imp)
  
summary(master_imputed) #correct


# Generate other variables 

master_imputed$imd_quintile <- ifelse(master_imputed$imd <= 4, 1,
                           ifelse(master_imputed$imd <= 8, 2,
                                  ifelse(master_imputed$imd <= 12, 3,
                                         ifelse(master_imputed$imd <= 16, 4,
                                                5))))
master_imputed <- master_imputed %>%
  mutate(imd_quintile=recode(imd_quintile,
                             `1`="Quintile 1(least deprived)",
                             `2`="Quintile 2",
                             `3`="Quintile 3",
                             `4`="Quintile 4",
                             `5`="Quintile 5"))

table(master_imputed$imd_quintile)


#  Generate age_group_uk variable according to ONS data

# Define a function to assign age groups
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

master_imputed <- master_imputed %>%
  mutate(age_group_uk=sapply(age, assign_age_group_uk))

master_imputed <- master_imputed %>%
  select(patid, imd, imd_quintile, gender, ethnic_2, age, age_group_uk, everything())

saveRDS(master_imputed, "D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds")


################################################################################
# 2. Generate infections datasets for each disease, 19 diseases                #
# (one row per consultation/visit)                                             #
################################################################################

# (1) Asthma & otitis media

# For asthma and otitis media, use Nam's new data

visit1621 <- fread("D:\\CPRDData\\Analysis_Nam\\Final datasets 19Jan\\0-visit-16-21-1May.txt")

visit19 <- visit1621 %>%
  mutate(obsdate = as.Date(obsdate)) %>%
  filter(obsdate >= as.Date("2019-01-01") & obsdate <= as.Date("2019-12-31"))

#----------------------------------------------
# Otitis media

visit19_media <- visit19 %>%
  filter(inf_otitis==1) %>%
  dplyr:: select(patid, pracid, obsdate, age, gender, abyes, ab_substance) #n=196583

visit19_media <- visit19_media %>%
  rename(date = obsdate)

saveRDS(visit19_media, "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/media.rds")

#----------------------------------------------
# Asthma 

visit19_asthma <- visit19 %>%
  filter(inf_asthma==1) %>%
  dplyr:: select(patid, pracid, obsdate, age, gender, abyes, ab_substance) #n=71370

visit19_asthma <- visit19_asthma %>%
  rename(date = obsdate)

saveRDS(visit19_asthma, "D:/CPRDData/Analysis_Murong/Data/Disease_imputed/asthma.rds")

# (2) Other diseases, use
infection <- readRDS("D:/CPRDData/Analysis_Murong/Data/infection_1621_full.rds") 

table(infection$inf_diagnosis)
table(infection$inf_position)

summary(infection$date) #2019-01-01 to 2019-12-31


# Rename abu to abyes to be consistent with Nam's new dataset visit1621

infection <- infection %>%
  rename(abyes = abu)


#----------------------------------------------
# COPD

copd <- infection %>%
  filter(inf_diagnosis == "COPD exacerbation")

#----------------------------------------------
# Lower RTI
lrti <- infection %>%
  filter(inf_diagnosis =="bronchinits " |
           inf_diagnosis == "lower RTI" | inf_diagnosis == "pneumoniae")

#----------------------------------------------
# Cough
# excluding COPD & asthma exacerbation, pneumoniae, and lower RTI if these diagnosis happened on the same visit for diagnosis with cough

cough <- infection %>%
  filter(inf_diagnosis %in% c("cough", "COPD exacerbation", "asthma exacerbation", "pneumoniae", "lower RTI"))

# Generate an indicator whether the visit for cough should be excluded
cough <- cough %>%
  group_by(patid, date) %>%
  mutate(exclude=if(
    any(inf_diagnosis %in% c("COPD exacerbation", "asthma exacerbation", "pneumoniae", "lower RTI") ) &
    any(inf_diagnosis == "cough")
  ) 1 else 0
  ) %>%
  ungroup()

table(cough$exclude)


check <- cough %>% filter(exclude==1) %>% dplyr::select(patid, date, inf_diagnosis, exclude)         
rm(check)

# Generate the cough dataset    
cough <- cough %>%
  filter(inf_diagnosis == "cough") %>%
  filter(exclude==0)

table(cough$inf_diagnosis)
table(cough$exclude)

#----------------------------------------------
# Upper RTI
urti <- infection %>%
  filter(inf_diagnosis == "sore throat" | inf_diagnosis == "upper RTI" | inf_diagnosis == "sinusitis")

#----------------------------------------------
# Individual respiratory infections

sore <- infection %>%
  filter(inf_diagnosis == "sore throat")

sinu <- infection %>%
  filter(inf_diagnosis == "sinusitis")

uurti <- infection %>%
  filter(inf_diagnosis == "upper RTI")

bron <- infection %>%
  filter(inf_diagnosis == "bronchinits")

ulrti <- infection %>%
  filter(inf_diagnosis == "lower RTI")

pneu <- infection %>%
  filter(inf_diagnosis == "pneumoniae")

#----------------------------------------------
# Gastroenteritis
gast <- infection %>%
  filter(inf_organ=="gastroenteritis") 

#----------------------------------------------
# Cellulite, impetigo, otitis externa

cell <- infection %>%
  filter(inf_diagnosis=="cellulitis")

impe <- infection %>%
  filter(inf_diagnosis=="impetigo")

exter <- infection %>%
  filter(inf_diagnosis=="otitis external")


#----------------------------------------------
# UTI & pyelonethritis

utio <- infection %>%
  filter(inf_diagnosis=="UTI") %>%
  filter(gender==2)

pyel <- infection %>%
  filter(inf_diagnosis=="pyelonethritis")%>%
  filter(gender==2)


#----------------------------------------------
# Prostatitis 
pros <- infection %>%
  filter(inf_diagnosis=="prostatitis")%>%
  filter(gender==1)


# Save dataset

setwd("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/")    

saveRDS(copd, "copd.rds") 
saveRDS(lrti, "lrti.rds") 
saveRDS(cough, "cough.rds") 
saveRDS(urti, "urti.rds")

saveRDS(gast, "gast.rds")
saveRDS(cell, "cell.rds")
saveRDS(impe, "impe.rds") 
saveRDS(exter, "exter.rds") 
saveRDS(utio, "utio.rds") 
saveRDS(pyel, "pyel.rds") 
saveRDS(pros, "pros.rds") 

saveRDS(sore, "sore.rds") 
saveRDS(sinu, "sinu.rds") 
saveRDS(uurti, "uurti.rds") 
saveRDS(bron, "bron.rds") 
saveRDS(ulrti, "ulrti.rds") 
saveRDS(pneu, "pneu.rds") 



################################################################################
# 3. Generate disease master dataset, number of new case & new abx             #
################################################################################

# Function 

generate_outcome <- function(infection){
  
  # Load master dataset
  master <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds")
  assign("master", master, envir = .GlobalEnv)
  
  
  # Construct file paths based on the infection name
  infection_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, ".rds") 
  
  # Open datasets
  infection_data <- readRDS(infection_path)
    
  # Convert the date column to Date format & arrange dataset by patid and consultation date
  infection_data <- infection_data %>%   
    mutate(date = as.Date(date)) %>%
    arrange(patid, date)
  
  #-----------------------------------------------------------------------------
  
  # (1) Generate number of new cases

  infection_new_case <- infection_data %>%   
    group_by(patid) %>%
    mutate(new_case = if_else(date - lag(date) > 30, 1, 0, missing = 0)) %>%  # generate variable
    ungroup() %>%
    group_by(patid) %>%
    mutate(new_case = if_else(row_number() == 1, 1, new_case)) %>%  
    summarise(n_new_case = sum(new_case)) %>%  
    ungroup()
  
  
  # Merge df dataset with master dataset to generate a master dataset for each disease
  infection_master <- master %>% 
    left_join(infection_new_case %>% dplyr::select(patid,n_new_case),by="patid") 
  
  #Replace NA with 0 for n_new_case
  infection_master <- infection_master %>%
    mutate(n_new_case = ifelse(is.na(n_new_case),0,n_new_case))

  
  #-----------------------------------------------------------------------------
  
  
  # (2) Generate number of new antibiotic prescription 
  
  infection_abx <- infection_data %>%
    filter(abyes==1) %>%      #includes only those with antibiotic prescription 
    group_by(patid) %>%
    mutate(new_abx = if_else(date - lag(date) > 30, 1, 0, missing = 0)) %>%  # generate variable
    ungroup() %>%
    group_by(patid) %>%
    mutate(new_abx = if_else(row_number() == 1, 1, new_abx)) %>%  
    summarise(n_new_abx = sum(new_abx)) %>% 
    ungroup()
  
  # Merge number of new abx  with disease master dataset
  infection_master <- infection_master %>% 
    left_join(infection_abx %>% dplyr::select(patid,n_new_abx),by="patid")
  
  infection_master <- infection_master %>% 
    mutate(n_new_abx = ifelse(is.na(n_new_abx),0,n_new_abx))
  
  #-----------------------------------------------------------------------------
  
  # (2) Generate master dataset for those who have at least one visit 
  
  infection_one_visit <- infection_data %>%
    group_by(patid) %>%
    slice_sample(n = 1) %>% #random select visit per patient
    ungroup()
  
  # and then merge with disease master dataset 
  
  infection_master <- infection_master  %>%
    left_join(infection_one_visit %>% dplyr::select(patid, abyes), by="patid") 
  
  infection_master_visit <- infection_master  %>%
    filter(n_new_case > 0)     # only includes patients with visit
  
  #-----------------------------------------------------------------------------
  
  # (3) Save datasets
  master_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master.rds") 
  master_visit_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_visit.rds") 

  saveRDS(infection_master, master_path)
  saveRDS(infection_master_visit, master_visit_path)
}




# Set parallel

infections <- c("asthma", "copd", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti")  # exclude "utio", "pyel", "pros" here as they are one-sex disease

num_cores <- 20

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterExport(cl, c("infections", "generate_outcome"))

results <- parLapply(cl, infections, generate_outcome)

beep("ping")


stopCluster(cl)

beep("ping")



# UTI & pyel

# Function 

generate_outcome_female <- function(infection){
  
  # Load master dataset
  master <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds")
  assign("master", master, envir = .GlobalEnv)
  
  # Only keep female for UTI 
  master <- master %>%
    filter(gender == "Female")
  
  # Construct file paths based on the infection name
  infection_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, ".rds") 
  
  # Open datasets
  infection_data <- readRDS(infection_path)
  
  # Convert the date column to Date format & arrange dataset by patid and consultation date
  infection_data <- infection_data %>%   
    mutate(date = as.Date(date)) %>%
    arrange(patid, date)
  
  #-----------------------------------------------------------------------------
  
  # (1) Generate number of new cases
  
  infection_new_case <- infection_data %>%   
    group_by(patid) %>%
    mutate(new_case = if_else(date - lag(date) > 30, 1, 0, missing = 0)) %>%  # generate variable
    ungroup() %>%
    group_by(patid) %>%
    mutate(new_case = if_else(row_number() == 1, 1, new_case)) %>%  
    summarise(n_new_case = sum(new_case)) %>%  
    ungroup()
  
  
  # Merge df dataset with master dataset to generate a master dataset for each disease
  infection_master <- master %>% 
    left_join(infection_new_case %>% dplyr::select(patid,n_new_case),by="patid") 
  
  #Replace NA with 0 for n_new_case
  infection_master <- infection_master %>%
    mutate(n_new_case = ifelse(is.na(n_new_case),0,n_new_case))
  
  
  #-----------------------------------------------------------------------------
  
  
  # (2) Generate number of new antibiotic prescription 
  
  infection_abx <- infection_data %>%
    filter(abyes==1) %>%      #includes only those with antibiotic prescription 
    #mutate(date = as.Date(date)) %>%
    group_by(patid) %>%
    mutate(new_abx = if_else(date - lag(date) > 30, 1, 0, missing = 0)) %>%  # generate variable
    ungroup() %>%
    group_by(patid) %>%
    mutate(new_abx = if_else(row_number() == 1, 1, new_abx)) %>%  
    summarise(n_new_abx = sum(new_abx)) %>%  
    ungroup()
  
  # Merge number of new abx  with disease master dataset
  infection_master <- infection_master %>% 
    left_join(infection_abx %>% dplyr::select(patid,n_new_abx),by="patid")
  
  infection_master <- infection_master %>% 
    mutate(n_new_abx = ifelse(is.na(n_new_abx),0,n_new_abx))
  
  #-----------------------------------------------------------------------------
  
  # (2) Generate master dataset for those who have at least one visit 
  
  infection_one_visit <- infection_data %>%
    group_by(patid) %>%
    slice_sample(n = 1) %>% 
    ungroup()
  
  # and then merge with disease master dataset 
  
  infection_master <- infection_master  %>%
    left_join(infection_one_visit %>% dplyr::select(patid, abyes), by="patid") 

  
  infection_master_visit <- infection_master  %>%
    filter(n_new_case > 0)     
  
  #-----------------------------------------------------------------------------
  
  # (3) Save datasets
  master_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master.rds") 
  master_visit_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_visit.rds") 
  
  saveRDS(infection_master, master_path)
  saveRDS(infection_master_visit, master_visit_path)
}




# Set parallel

infections <- c("utio", "pyel")


num_cores <- 5

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterExport(cl, c("infections", "generate_outcome_female"))

results <- parLapply(cl, infections, generate_outcome_female)

beep("ping")


stopCluster(cl)

beep("ping")




# Pros

# Function 

generate_outcome_male <- function(infection){
  
  # Load master dataset
  master <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds")
  assign("master", master, envir = .GlobalEnv)
  
  # Only keep female for UTI 
  master <- master %>%
    filter(gender == "Male")
  
  # Construct file paths based on the infection name
  infection_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, ".rds") 
  
  # Open datasets
  infection_data <- readRDS(infection_path)
  
  # Convert the date column to Date format & arrange dataset by patid and consultation date
  infection_data <- infection_data %>%   
    mutate(date = as.Date(date)) %>%
    arrange(patid, date)
  
  #-----------------------------------------------------------------------------
  
  # (1) Generate number of new cases
  
  infection_new_case <- infection_data %>%   
    group_by(patid) %>%
    mutate(new_case = if_else(date - lag(date) > 30, 1, 0, missing = 0)) %>%  # generate variable
    ungroup() %>%
    group_by(patid) %>%
    mutate(new_case = if_else(row_number() == 1, 1, new_case)) %>%  
    summarise(n_new_case = sum(new_case)) %>%  
    ungroup()
  
  
  # Merge df dataset with master dataset to generate a master dataset for each disease
  infection_master <- master %>% 
    left_join(infection_new_case %>% dplyr::select(patid,n_new_case),by="patid") 
  
  #Replace NA with 0 for n_new_case
  infection_master <- infection_master %>%
    mutate(n_new_case = ifelse(is.na(n_new_case),0,n_new_case))
  
  
  #-----------------------------------------------------------------------------
  
  
  # (2) Generate number of new antibiotic prescription 
  
  infection_abx <- infection_data %>%
    filter(abyes==1) %>%      #includes only those with antibiotic prescription 
    #mutate(date = as.Date(date)) %>%
    group_by(patid) %>%
    mutate(new_abx = if_else(date - lag(date) > 30, 1, 0, missing = 0)) %>%  
    ungroup() %>%
    group_by(patid) %>%
    mutate(new_abx = if_else(row_number() == 1, 1, new_abx)) %>% 
    summarise(n_new_abx = sum(new_abx)) %>% 
    ungroup()
  
  # Merge number of new abx  with disease master dataset
  infection_master <- infection_master %>% 
    left_join(infection_abx %>% dplyr::select(patid,n_new_abx),by="patid")
  
  infection_master <- infection_master %>% 
    mutate(n_new_abx = ifelse(is.na(n_new_abx),0,n_new_abx))
  
  #-----------------------------------------------------------------------------
  
  # (2) Generate master dataset for those who have at least one visit 
  
  infection_one_visit <- infection_data %>%
    group_by(patid) %>%
    slice_sample(n = 1) %>% 
    ungroup()
  
  # and then merge with disease master dataset 
  
  infection_master <- infection_master  %>%
    left_join(infection_one_visit %>% dplyr::select(patid, abyes), by="patid") 
  
  infection_master_visit <- infection_master  %>%
    filter(n_new_case > 0)    
  
  #-----------------------------------------------------------------------------
  
  # (3) Save datasets
  master_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master.rds") 
  master_visit_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_visit.rds") 
  
  saveRDS(infection_master, master_path)
  saveRDS(infection_master_visit, master_visit_path)
}



# Set parallel

infections <- c("pros")

num_cores <- 5

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterExport(cl, c("infections", "generate_outcome_male"))

results <- parLapply(cl, infections, generate_outcome_male)

beep("ping")

stopCluster(cl)

beep("ping")




################################################################################
# 4. Generate master datasets for children and adults                          #
################################################################################

update_sample <- function(infection) {
  
  # ------     Generate master dataset & dataset for children/adults      -----#
  
  # Master dataset 
  master_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master.rds") 
  
  infection_master <- readRDS(master_path)
  
  # Generate master datasets for children and adults 
  
  df_child <- infection_master %>%
    filter(age < 18)
  
  df_adult <- infection_master %>%
    filter(age >= 18)
  
  # Construct file names
  
    file_name_child <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_child.rds")
    file_name_adult <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_adult.rds")
  
  # Save datasets  
  saveRDS(df_child, file_name_child)
  saveRDS(df_adult, file_name_adult)
  
  
  
  # --------    Generate master_visit dataset & by children/adults     --------#
  
  # Master dataset 
  master_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_visit.rds") 
  
  infection_master <- readRDS(master_path)
  
  # Generate master datasets for children and adults 
  
  df_child <- infection_master %>%
    filter(age < 18)
  
  df_adult <- infection_master %>%
    filter(age >= 18)
  
  # Construct file names
  
  file_name_child <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_visit_child.rds")
  file_name_adult <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, "_master_visit_adult.rds")
  
  # Save datasets  
  saveRDS(df_child, file_name_child)
  saveRDS(df_adult, file_name_adult)
  
}


infections_list <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "media", "exter", "gast", "utio", "pyel", "pros",
                     "sore", "sinu", "uurti", "bron", "ulrti")  


num_cores <- detectCores()-1

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterExport(cl, c("infections_list", "update_sample"))

results <- parLapply(cl, infections_list, update_sample)

stopCluster(cl)

beep("ping")


















