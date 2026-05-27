# GAMLSS for adult/child dataset

# Load necessary libraries
library(dplyr)
library(gamlss)
library(parallel)
library(doParallel)
library(beepr)
library(profvis)

################################################################################
# 1. Run parallel                                                              #
################################################################################

#------------------------------------------------------------------------------#
# CHILDREN                                                                     #
#------------------------------------------------------------------------------#

# (1) Create a function 

# Update gamlss_abx_child_uti_resp function with tryCatch to handle errors

gamlss_abx_child_uti <- function(infection, sample_size) {
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
      dplyr::select(patid, bs_age_1, bs_age_2, bs_age_3, bs_age_4, gender, ethnic_2, n_new_abx, bs_imd_1, bs_imd_2, log_timeatrisk_abx, vacc_pneu, vacc_flu)
    
    infection_sample <- infection_sample %>%
      filter(complete.cases(.))
    
    sample_time <- Sys.time()  
    logs <- c(logs, paste0("Sample selection finish: ", sample_time))
    
    
    start_model <- Sys.time()
    logs <- c(logs, paste0("GAMLSS start: ", start_model))
    
    # # Run GAMLSS models
    
    m1_infection <- gamlss(n_new_abx ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + ethnic_2 +
                             offset(log_timeatrisk_abx),
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 2/Model/Child/", "m1_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!! m0/1/...
    saveRDS(m1_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!! 
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))  
    })
}

# (2) Set up parallel

# Infections for child
infections <- c("utio", "pyel")

num_cores <- 5

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_abx_child_uti", "infections"))

# (3) Run GAMLSS  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- parLapply(cl, infections, function(infection) gamlss_abx_child_uti(infection, sample_size = 1))

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")

print(results)

#------------------------------------------------------------------------------#
# ADULTS                                                                       #
#------------------------------------------------------------------------------#

# (1) Function 


gamlss_abx_adult_one_sex <- function(infection, sample_size) {
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
      dplyr::select(patid, bs_age_1, bs_age_2, bs_age_3, bs_age_4, gender, ethnic_2, n_new_abx, bs_imd_1, bs_imd_2, log_timeatrisk_abx, vacc_pneu, vacc_flu)
    
    infection_sample <- infection_sample %>%
      filter(complete.cases(.))
    
    sample_time <- Sys.time()  
    logs <- c(logs, paste0("Sample selection finish: ", sample_time))
    
    
    start_model <- Sys.time()
    logs <- c(logs, paste0("GAMLSS start: ", start_model))
    
    # # Run GAMLSS models
    
    m1_infection <- gamlss(n_new_abx ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + ethnic_2 +
                             offset(log_timeatrisk_abx),
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))

    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 2/Model/Adult/", "m1_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m1_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))  
    })
}



# (2) Set up parallel

# Infections for ADULTS

infections <- c("utio", "pyel", "pros")

num_cores <- 10

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_abx_adult_one_sex", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_abx_adult_one_sex(infection, sample_size = 1))

stopCluster(cl)

beep("ping")

print(results)








