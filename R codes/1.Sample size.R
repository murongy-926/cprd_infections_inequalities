# Sample selection flow chart  

setwd("D:/CPRDData/Analysis_Murong/Data/")    
dir()

install.packages("data.table")
install.packages("tidyverse")
install.packages("dplyr")

library(bit64)
library(data.table)
library(tidyverse)
library(dplyr)

################################################################################
# 0. Study population extracted in our project                                 #
################################################################################

pat_visit_1619 <- readRDS("D:/CPRDData/Analysis_Murong/Data/pat_visit_1619.rds")


################################################################################
# 1.Generate a master dataset  (one row per patient)                           #
################################################################################

# 2.1 Data preparation 

# Check practice last collection date

practice <- fread("D:/CPRDData/data1/202309_CPRDAurum/202309_CPRDAurum_Practices.txt")
practice %>%
  count(pracid)

practice$lcd <- as.Date(practice$lcd,format="%d/%m/%Y")

summary(practice$lcd) # Min=2018-02-22

# Generate an indicator

practice$present <- ifelse(practice$lcd<as.Date("2019-01-01"),0,1)
table(practice$present)  

# 2.2 Generate master dataset

practice_selected <- practice %>% dplyr::select(pracid,lcd,present)
master <- left_join(pat_visit_1619, practice_selected, by="pracid") 

table(master$present) 
summary(master$pracid)

rm(practice_selected)


# Collect patient registration start and end date 

patient <-fread("D:/CPRDData/data1/202309_CPRDAurum/202309_CPRDAurum_AcceptablePats.txt")    

patient_selected <- patient %>% dplyr::select(-c("uts","lcd","region"))

master <- inner_join(master,patient_selected,by=c("patid","pracid")) 

rm(patient,patient_selected)

master <- master %>%
  rename(prac_present=present)

master$regstartdate <- as.Date(master$regstartdate,format="%d/%m/%Y")

master$regenddate <- as.Date(master$regenddate,format="%d/%m/%Y")

# Generate an indicator

master$pat_present <- ifelse(master$regenddate<as.Date("2019-01-01"),0,1)
table(master$pat_present)    
summary(master$pat_present) 

master_prac_0 <- master %>%
  filter(prac_present==0)
table(master_prac_0$pat_present)

rm(master_prac_0)

# Copy master dataset
master_full <- master
saveRDS(master_full, file="master_full_1619_update.rds")

master<- readRDS("master_full_1619_update.rds")




################################################################################
# 1. Validated CPRD record - Exclusion criteria
################################################################################

# 1.1 Validate data according to CPRD and NDPH guidance

# Exclude 29 practices that appear likely to have merged into other contributing practices, that CPRD suggests to exclude from the study
gp.exc <- data.frame(pracid = as.integer(c(20024, 20036, 20091, 20202, 20254, 20389, 20430, 20469, 20487, 20552, 20554, 20734, 20790, 20803, 20868, 20996, 21001, 21078, 21118, 21172, 21173, 21277, 21334, 21390, 21444, 21451, 21553, 21558, 21585)))
master <- master %>% anti_join(gp.exc, by="pracid") # Drop = 10443127-10364831=78296      

# Examine patient category and exclude non-permanent patients if exist 
summary(master$patienttypeid)

# Exclude visits with gender not identified as male or female

master <- master %>%
  filter(gender != 3) 

# (5) Exclude visits with empty registration date 
summary(master$regstartdate) 

# Exclude visits with registration date < 01/01/1900  #Drop = 0

# Exclude visits with a date before the year of birth (done in infections dataset)

# Exclude visits with a date later than the date of death (done in infections dataset)

# Exclude visits with a date later than the practice last collection date (done in infections dataset) 

# Exclude visits with the registration start date later than the practice last collection date

master <- master %>%
  filter(regstartdate < lcd) 

# (11) Exclude visits with empty year of birth 
summary(is.na(master$yob)) #0 missing


# (1) Exclude the visit corresponding to a patient aged more than 115 years old #Drop = 0 

current_year <- 2019 #Consultation year is 2019

master <- master %>%
  mutate(age=current_year-yob)

summary(master$age)

master <- master %>%
  filter(age <=115) 


################################################################################
#  In this project with research period in 2019: exclude patients/practices that are not registered on/after 2019-01-01
################################################################################

#(1).practice present at least on/after 2019-01-01,prac_present==1
master <- master %>%
  filter(!(prac_present==0)) 

#(2).patient registers at least on/after 2019-01-01,pat_present==0
master <- master %>%
  filter(pat_present ==1 | is.na(pat_present)) 

#(3) Remove patients born after 2019.12.31
summary(master$yob)

rows_to_delete <- which(master$yob>2019)
master <- master[-rows_to_delete]  

master <- master %>%
  filter(yob<2020) 

# (4) Remove patients registered after 2019.12.31
# Details see Section 3.1 below

# 3. Save dataset

saveRDS(master, file="master.rds")




################################################################################
# 3. Generate covariates in master dataset                                      #
################################################################################

# 3.1 Generate registration start and end date ####

# 1. Registration start date

first_date <- readRDS("First date/first date.rds") #first observation date for patients in CPRD 

summary(is.na(first_date$obsdate)) #No missing

master <- master %>%
  left_join(first_date %>% select(-pracid) %>% rename(first_obsdate = obsdate), by="patid")

summary(master$first_obsdate)

summary(master$regstartdate)

master <- master %>%
  mutate(First_obsdate = as.Date(first_obsdate),  # Ensure first_obsdate is a Date object
         First_obsdate = dplyr::if_else(first_obsdate < as.Date("1987-01-01"), 
                                        as.Date("1987-01-01"), first_obsdate)) 
summary(master$First_obsdate)

master <- master %>%
  mutate(Regstartdate = as.Date(regstartdate),  
         Regstartdate = dplyr::if_else(regstartdate < as.Date("1987-01-01"), 
                                       as.Date("1987-01-01"), regstartdate)) 

summary(master$Regstartdate)

# Generate registration start date
# Rule: min betwen regstartdate & first_obsdate

master <- master %>%
  mutate(startdate=pmin(Regstartdate, First_obsdate, na.rm=TRUE)) # Set na.rm=TRUE to ignore missing values in First_obsdate, pmin() will choose a valid value, i.e. use Regstartdate 

summary(master$startdate)

# Exclude patients registered after 2019.12.31
master <- master %>%
  filter(startdate<as.Date("2020-01-01"))  

# 2. Registration end date
# Rule: minimum between practice last collection date and patient registration end date

master <- master %>%
  mutate(enddate=pmin(lcd,regenddate,na.rm=TRUE))  
master_2019 <- master %>%
  filter(lcd<as.Date("2019-12-31") & lcd>as.Date("2019-01-01"))  
rm(master_2019)



# 3.2 Indicator of Whole year(=2019) present ####

master$wholeyear_present <- ifelse((master$startdate<=as.Date("2019-01-01") & master$enddate>=as.Date("2019-12-31")),1,0)
summary(master$wholeyear_present)  
table(master$wholeyear_present)   




# 3.3 IMD 

# 1.1 IMD data ####

imd <- fread("D:/CPRDData/23_003072_Aurum_1/patient_2019_imd_23_003072.txt") 

master <- master %>%
  left_join(imd %>% select(-pracid), by = "patid")

summary(is.na(master$e2019_imd_20))


# 3.4 Age, gender and ethnicity ####

# 1. Generate age_group variable by patient age

# Generate age variable 

current_year <- 2019 #Consultation year is 2019

master <- master %>%
  mutate(age=current_year-yob)


master$imd_quintile <- ifelse(master$e2019_imd_20 <= 4, 1,
                              ifelse(master$e2019_imd_20 <= 8, 2,
                                     ifelse(master$e2019_imd_20 <= 12, 3,
                                            ifelse(master$e2019_imd_20 <= 16, 4,
                                                   5))))

master <-master %>%
  select(patid,pracid,gender,age,yob,mob,e2019_imd_20,imd_quintile,everything())

# Define a function to assign age groups
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

master <- master %>%
  mutate(age_group=sapply(age, assign_age_group))

# 2. Gender
master$gender <- as.factor(master$gender)

master$gender <- factor(master$gender,
                        levels = c(1, 2),
                        labels = c("Male", "Female"))

master$gender <- relevel(factor(master$gender), ref = "Male")

# 3. Ethnicity 
ethnic <- readRDS("ethnic.rds")

master <- master %>%
  left_join(ethnic %>% dplyr::select(patid, ethnic_2, ethnic_4), by = "patid")  

summary(is.na(master$ethnic_4))
summary(is.na(master$ethnic_2))

#Reorder the columns in the master dataset
master <-master %>%
  rename(imd = e2019_imd_20) %>%
  dplyr::select(patid,pracid,gender,age,age_group,yob,mob,ethnic_2, ethnic_4, imd,imd_quintile,everything())

# Save dataset 
saveRDS(master,"master_1619_with_missing_covar.rds") 



################################################################################
# 4.Merge contributors into master and disease master datasets                  #
################################################################################

# Covariates: 
# smoking (has already merged to master dataset)
# hazardous & harmful alcohol drinking
# BMI
# vaccination - pneumococcal, flu, rotavirus


# 6.1 Open covariates datasets ###
smk_single <- readRDS("Covariates/Smoking/smoking_single_update.rds")
alcohol <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/Alcohol/alcohol_recent.rds")
bmi <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/BMI/bmi.rds")
pneu <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/pneumococcal/pneu_sample.rds")
flu <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/flu/flu_sample.rds")
rotav <- readRDS("D:/CPRDData/Analysis_Murong/Data/Covariates/Vaccination/rotavirus/rotav_sample.rds")

# 6.2  Merge into master dataset ###

# (0) Smoking

master <- master %>%
  left_join(smk_single %>% select(patid, current_smkstatus_new), by="patid")

sum(is.na(master$current_smkstatus_new)) 

table(master$current_smkstatus_new) 

# Rename variable
master <- master %>%
  rename(current_smkstatus=current_smkstatus_new)

# Replace NA as non-smoker

master <- master %>%
  mutate(current_smkstatus_update=ifelse(is.na(current_smkstatus),"non-smoker",current_smkstatus))

# (1) alcohol

summary(is.na(alcohol)) 

master <- master %>%
  left_join(alcohol %>%
              rename(term_alcohol = Term) %>%
              dplyr::select(patid, term_alcohol),
            by = "patid"
  )


# Generate variable

master <- master %>%
  mutate(alcohol = ifelse(!is.na(term_alcohol), 1, 0)) #1 indicates that a patient is in the alcohol dataset and 0 indicates otherwise

table(master$alcohol) 

# (2) BMI

master <- master %>%
  left_join(bmi %>% dplyr::select(patid, bmi),
            by = "patid")

master$bmi <- round(master$bmi, 2)


summary(is.na(master$bmi))


# (3) Vaccination 

# pneucomoccal ###

# Merge

summary(is.na(pneu))

pneu <- pneu %>%
  mutate(pneu_indicator = 1) 


master <- master %>%
  left_join(pneu %>%
              dplyr::select(patid, pneu_indicator),
            by = "patid"
  )

# Generate variable

master <- master %>%
  mutate(vacc_pneu = ifelse(!is.na(pneu_indicator), 1, 0)) 

table(master$vacc_pneu) # 


# flu ###

# Merge 

summary(is.na(flu)) 

flu <- flu %>%
  mutate(flu_indicator = 1) 

master <- master %>%
  left_join(flu %>%
              dplyr::select(patid, flu_indicator),
            by = "patid"
  )

# Generate variable

master <- master %>%
  mutate(vacc_flu = ifelse(!is.na(flu_indicator), 1, 0)) 

table(master$vacc_flu) 


saveRDS(master,"master_1619_with_missing_covar.rds") #N=9426835


################################################################################
# 7 Drop missing covariates
################################################################################

master_nomiss <- master %>%
  filter(!is.na(imd))

master_nomiss <- master_nomiss %>%
  filter(!is.na(ethnic_2))


master_nomiss <- master_nomiss %>%
  filter(current_smkstatus_update != "M") 

master_nomiss <- master_nomiss %>%
  filter(!(is.na(bmi) & age >= 18)) 
















