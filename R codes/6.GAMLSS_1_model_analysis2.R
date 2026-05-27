# GAMLSS for adult/child dataset
# Analysis 2: inequality in antibiotic prescription 

# Load necessary libraries
library(dplyr)
library(gamlss)
library(parallel)
library(doParallel)
library(beepr)
library(profvis)
library(bit64)

################################################################################
# 0. Check if variables are available / in right class in the datasets         #
################################################################################

# Already checked in analysis 1

################################################################################
# 2. Run parallel                                                              #
################################################################################


# 2.1 Parallel by diseases                                                     

#------------------------------------------------------------------------------#
# CHILDREN                                                                     #
#------------------------------------------------------------------------------#

# (1) Create a function 

# Update gamlss_abx_child_resp function with tryCatch to handle errors

gamlss_abx_child <- function(infection, sample_size) {
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
      dplyr::select(patid, bs_age_1, bs_age_2, bs_age_3, bs_age_4, gender, ethnic_2, bs_imd_1, bs_imd_2, vacc_pneu, vacc_flu,
                    n_new_abx, log_timeatrisk_abx)
    
    infection_sample <- infection_sample %>%
      filter(complete.cases(.))
    
    sample_time <- Sys.time()  
    logs <- c(logs, paste0("Sample selection finish: ", sample_time))
    
    
    start_model <- Sys.time()
    logs <- c(logs, paste0("GAMLSS start: ", start_model))
    
    # # Run GAMLSS models

    m1_infection <- gamlss(n_new_abx ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk_abx),
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
                               
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 2/Model/Child/", "m1_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!! m0/1/...
    saveRDS(m1_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!! m0/1/...
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    
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
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast",
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude uti and pyel, do it seperately as UTI only includes female

num_cores <- 5

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_abx_child", "infections"))

# (3) Run GAMLSS  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#results <- parLapply(cl, infections, gamlss_abx_child)
results <- parLapply(cl, infections, function(infection) gamlss_abx_child(infection, sample_size = 1))

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")

print(results)

#------------------------------------------------------------------------------#
# ADULTS                                                                       #
#------------------------------------------------------------------------------#

# (1) Function 


gamlss_abx_adult <- function(infection, sample_size) {
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
      dplyr::select(patid, bs_age_1, bs_age_2, bs_age_3, bs_age_4, gender, ethnic_2,bs_imd_1, bs_imd_2, vacc_pneu, vacc_flu, 
                    n_new_abx, log_timeatrisk_abx)
    
    infection_sample <- infection_sample %>%
      filter(complete.cases(.))
    
    sample_time <- Sys.time()  
    logs <- c(logs, paste0("Sample selection finish: ", sample_time))
    
    
    start_model <- Sys.time()
    logs <- c(logs, paste0("GAMLSS start: ", start_model))
    
    # # Run GAMLSS models

    m1_infection <- gamlss(n_new_abx ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk_abx),
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
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

infections <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast",
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude "utio", "pyel", "pros", do them separately as they only includes one sex group


num_cores <- 5

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_abx_adult", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_abx_adult(infection, sample_size = 1))

stopCluster(cl)

beep("ping")

print(results)







