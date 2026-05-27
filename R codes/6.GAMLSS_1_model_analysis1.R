# GAMLSS for adult/child dataset
# Use model specification 1 only in this script (i.e. control for age sex ethnicity)

# Load necessary libraries
library(dplyr)
library(gamlss)
library(parallel)
library(doParallel)
library(beepr)
library(profvis)

################################################################################
# 0. Check if variables are available / in right class in the datasets         #
################################################################################

#------------------------------------------------------------------------------#
# (1) Availability of variables                                                #
#------------------------------------------------------------------------------#

library(foreach)
library(doParallel)

# Define all infections
infections <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "media", "exter", "gast", "utio", "pyel", "pros",
                "sore", "sinu", "uurti", "bron", "ulrti")

# Infections for CHILDREN
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", "utio", 
                "pyel", "sore", "sinu", "uurti", "bron", "ulrti")

# Infections for ADULTS

infections <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", "utio", "pyel", "pros",
                "sore", "sinu", "uurti", "bron", "ulrti")

# Define the path and the variables to check
master_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Data/Adult/"   #UPDATE Child/Adult accordingly

# variables for children #
required_vars <- c("patid", "age", "gender", "ethnic_2", 
                   "n_new_case", "n_new_abx","imd", "log_timeatrisk", "log_timeatrisk_abx", "vacc_pneu", "vacc_flu",  "abyes", 
                   "bs_imd_1", "bs_imd_2", 
                   "bs_age_1",  "bs_age_2",  "bs_age_3", "bs_age_4"  
                   ) 

# variables for adults #
required_vars <- c("patid", "age", "gender", "ethnic_2", 
                   "n_new_case", "n_new_abx","imd", "log_timeatrisk", "log_timeatrisk_abx", "vacc_pneu", "vacc_flu",  "abyes", 
                   "bs_imd_1", "bs_imd_2", 
                   "bs_age_1",  "bs_age_2",  "bs_age_3", "bs_age_4",
                   "bs_bmi_1", "bs_bmi_2", "bs_bmi_3","bs_bmi_4","bs_bmi_5",
                   "current_smkstatus_update", "bmi"
) 
  

# Function to load dataset and check variables
check_variables <- function(infection, master_path, required_vars) {
  file_path <- paste0(master_path, infection, "_master_adult.rds") ## UPDATE NAME ACCORDINGLY## _master_child/adult; _master_visit_child/adult
  if (!file.exists(file_path)) {
    return(paste("File not found for", infection))
  }
  
  dataset <- tryCatch(readRDS(file_path), 
                      error = function(e) return(paste("Error loading file for", infection)))
  
  if (is.character(dataset)) {
    return(dataset) # Return the error message if loading failed
  }
  
  missing_vars <- setdiff(required_vars, names(dataset))
  if (length(missing_vars) == 0) {
    return(paste("All required variables are present in", infection))
  } else {
    return(paste("Missing variables in", infection, ":", paste(missing_vars, collapse = ", ")))
  }
}

# Parallel setup
num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Check variables in parallel
results <- foreach(infection = infections, .combine = c, .packages = "base") %dopar% {
  check_variables(infection, master_path, required_vars)
}

# Stop the cluster
stopCluster(cl)

# Combine results with infection names
names(results) <- infections

# Print the results
print(results)

#------------------------------------------------------------------------------#
# (2) Check the class of variables                                             #
#------------------------------------------------------------------------------#

# Define the path and the variables to check
master_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Data/Adult/" ## UPDATE NAME ACCORDINGLY##
required_vars <- c("gender", "ethnic_2", "vacc_pneu", "vacc_flu", "current_smkstatus_update") #For adults
required_vars <- c("gender", "ethnic_2", "vacc_pneu", "vacc_flu") #For children

# Function to load dataset and check variables
check_variables <- function(infection, master_path, required_vars) {
  file_path <- paste0(master_path, infection, "_master_adult.rds") ## UPDATE NAME ACCORDINGLY##
  if (!file.exists(file_path)) {
    return(paste("File not found for", infection))
  }
  
  dataset <- tryCatch(readRDS(file_path), 
                      error = function(e) return(paste("Error loading file for", infection)))
  
  if (is.character(dataset)) {
    return(dataset) # Return the error message if loading failed
  }
  
  # Initialize a result variable to collect information about non-factor variables
  non_factor_vars <- c()
  
  # Loop through required_vars and check if they are factors
  for (var in required_vars) {
    if (!is.factor(dataset[[var]])) {
      non_factor_vars <- c(non_factor_vars, var)  # Add the variable name to non_factor_vars if not a factor
    }
  }
  
  # Return result based on whether any non-factor variables were found
  if (length(non_factor_vars) == 0) {
    return(paste("All required variables are factor variables in", infection))
  } else {
    return(paste("Non-factor variables in", infection, ":", paste(non_factor_vars, collapse = ", ")))
  }
}

# Parallel setup
num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Check variables in parallel
results <- foreach(infection = infections, .combine = c, .packages = "base") %dopar% {
  check_variables(infection, master_path, required_vars)
}

stopCluster(cl)

# Combine results with infection names
names(results) <- infections

# Print the results
print(results)

#~~~~~~~~~~~~#
#~~~~~~~~~~~~#

# Update variable to factor variable if necessary

# Go back to previous syntax and check if the variable is factor variable now 


################################################################################
# 1. Generate log_timeatrisk and unplinalised spline for IMD, Age and BMI      #
################################################################################

# Not necessary for total estimated group 


################################################################################
# 2. Run parallel                                                              #
################################################################################


# 2.1 Parallel by diseases                                                     

#------------------------------------------------------------------------------#
# CHILDREN                                                                     #
#------------------------------------------------------------------------------#

# (1) Create a function 

# Update gamlss_child_resp function with tryCatch to handle errors

gamlss_child <- function(infection, sample_size) {
  tryCatch({
    
    logs <- character()
    
    start_time <- Sys.time()
    logs <- c(logs, paste0("Start time: ", start_time))
    
    
    master_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Data/Child/", infection, "_master_child.rds")
    infection_master <- readRDS(master_path)
    
    read_time <- Sys.time() 
    logs <- c(logs, paste0("Read data finish: ", read_time))
    
    set.seed(123)  
    
    start_sample <- Sys.time()
    logs <- c(logs, paste0("Sample selection start: ", start_sample))
    
    infection_sample <- sample_frac(infection_master, sample_size) %>%
      dplyr::select(patid, bs_age_1, bs_age_2, bs_age_3, bs_age_4, gender, ethnic_2, n_new_case, bs_imd_1, bs_imd_2, log_timeatrisk, vacc_pneu, vacc_flu)
    
    infection_sample <- infection_sample %>%
      filter(complete.cases(.))
    
    sample_time <- Sys.time()  
    logs <- c(logs, paste0("Sample selection finish: ", sample_time))
    
    
    start_model <- Sys.time()
    logs <- c(logs, paste0("GAMLSS start: ", start_model))
    
    # # Run GAMLSS models
    
    m1_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk),
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))

    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Child/", "m1_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m1_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    return(list(error = paste("Error running GAMLSS for", infection, ": ", e$message)))
  })
}

# (2) Set up parallel

# Infections for child

infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude utio and pyel, do it seperately as UTI & pyel only includes female

num_cores <- 20

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_child", "infections"))

# (3) Run GAMLSS  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#results <- parLapply(cl, infections, gamlss_child)
results <- parLapply(cl, infections, function(infection) gamlss_child(infection, sample_size = 1))

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")



#------------------------------------------------------------------------------#
# ADULTS                                                                       #
#------------------------------------------------------------------------------#

# (1) Function


gamlss_adult <- function(infection, sample_size) {
  tryCatch({
    
    logs <- character()
    
    start_time <- Sys.time()
    logs <- c(logs, paste0("Start time: ", start_time))
    
    
    master_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Data/Adult/", infection, "_master_adult.rds")
    infection_master <- readRDS(master_path)
    
    read_time <- Sys.time() 
    logs <- c(logs, paste0("Read data finish: ", read_time))
    
    set.seed(123)  
    
    start_sample <- Sys.time()
    logs <- c(logs, paste0("Sample selection start: ", start_sample))
    
    infection_sample <- sample_frac(infection_master, sample_size) %>%
      dplyr::select(patid, bs_age_1, bs_age_2, bs_age_3, bs_age_4, gender, ethnic_2, n_new_case, bs_imd_1, bs_imd_2, log_timeatrisk, 
                    vacc_pneu, vacc_flu, current_smkstatus_update, bs_bmi_1, bs_bmi_2, bs_bmi_3, bs_bmi_4, bs_bmi_5)
    
    infection_sample <- infection_sample %>%
      filter(complete.cases(.))
    
    sample_time <- Sys.time()  
    logs <- c(logs, paste0("Sample selection finish: ", sample_time))
    
    
    start_model <- Sys.time()
    logs <- c(logs, paste0("GAMLSS start: ", start_model))
    
    # # Run GAMLSS models

    m1_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk),
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
                           
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Adult/", "m1_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m1_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    
                logs = logs))       
    
    
  },  error = function(e) {
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))
  })
}



# (2) Set up parallel

# Infections for ADULTS
infections <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude "utio", "pyel", "pros", do them separately as they only includes one sex group

num_cores <- 20

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_adult", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_adult(infection, sample_size = 1))

stopCluster(cl)

beep("ping")

print(results)




