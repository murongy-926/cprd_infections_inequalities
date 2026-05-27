# Analysis 1 model specification for models other than m1

################################################################################
################################################################################

# Load necessary libraries
library(dplyr)
library(gamlss)
library(parallel)
library(doParallel)
library(beepr)
library(profvis)

#------------------------------------------------------------------------------#
# CHILDREN
#------------------------------------------------------------------------------#

# Other model specifications, m0, m2-m3 


#---------------------------#
#---------------------------#

# m0 

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
    
    m0_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +offset(log_timeatrisk),
                           sigma.fo = ~ bs_imd_1 + bs_imd_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Child/", "m0_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m0_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
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
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude uti, do it seperately as UTI only includes female

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
print(results)

#---------------------------#
#---------------------------#

# m2

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
    
    m2_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk) + vacc_flu,
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Child/", "m2_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m2_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
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
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude uti, do it seperately as UTI only includes female

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
print(results)

#---------------------------#
#---------------------------#

# m3

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
    
    m3_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk) + vacc_pneu,
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Child/", "m3_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m3_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,     ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    return(list(error = paste("Error running GAMLSS for", infection, ": ", e$message)))
  })
}

# (2) Set up parallel

# Infections for child
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude uti, do it seperately as UTI only includes female

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
print(results)

beep("ping")


#------------------------------------------------------------------------------#
# ADULTS
#------------------------------------------------------------------------------#

# Other models m0, m2-m6

# m0 ###

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
    
    m0_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +offset(log_timeatrisk),
                           sigma.fo = ~ bs_imd_1 + bs_imd_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Adult/", "m0_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m0_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))  })
}



# (2) Set up parallel

# Infections for ADULTS
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude "utio", "pyel", "pros", do them separately as they only includes one sex group

# exclude copd as it is very slow, do it seperately


num_cores <- 25

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_adult", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_adult(infection, sample_size = 1))

stopCluster(cl)


#---------------------------#
#---------------------------#

# m2 ###


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
    
    m2_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk) + vacc_flu,
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Adult/", "m2_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m2_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))  })
}



# (2) Set up parallel

# Infections for ADULTS
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti")  #exclude "utio", "pyel", "pros", do them separately as they only includes one sex group

# exclude copd as it is very slow, do it seperately

num_cores <- 25

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_adult", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_adult(infection, sample_size = 1))

stopCluster(cl)

#---------------------------#
#---------------------------#

# m3 ###


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
    
    
    m3_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk) + vacc_pneu,
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Adult/", "m3_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m3_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))  })
}



# (2) Set up parallel

# Infections for ADULTS
infections <- c("asthma",  "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti")  #exclude "utio", "pyel", "pros", do them separately as they only includes one sex group

# exclude copd as it is very slow, do it seperately

num_cores <- 25

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_adult", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_adult(infection, sample_size = 1))

stopCluster(cl)

#---------------------------#
#---------------------------#

# m4 ###

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
    
    m4_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk) + current_smkstatus_update,
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Adult/", "m4_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m4_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
  
  
  },  error = function(e) {
    # Store the error message
    logs$error <- conditionMessage(e)
    return(list(success = FALSE, logs = logs))  })
}



# (2) Set up parallel

# Infections for ADULTS
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude "utio", "pyel", "pros", do them separately as they only includes one sex group

# exclude copd as it is very slow, do it seperately

num_cores <- 25

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_adult", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_adult(infection, sample_size = 1))

stopCluster(cl)

beep("ping")

#---------------------------#
#---------------------------#


#---------------------------#
#---------------------------#

# m6 ###

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
    m6_infection <- gamlss(n_new_case ~ bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2 +
                             offset(log_timeatrisk) + bs_bmi_1 + bs_bmi_2 + bs_bmi_3 + bs_bmi_4 + bs_bmi_5,
                           sigma.fo = ~  bs_imd_1 + bs_imd_2 +
                             bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + gender + ethnic_2,
                           family = NBI,
                           data = infection_sample,
                           method = RS(),
                           control = gamlss.control(c.crit = 0.001, n.cyc = 40))
    
    
    model_time <- Sys.time() 
    logs <- c(logs, paste0("GAMLSS end: ", model_time))
    
    # Save model                                                                
    
    save_model_path <- paste0("D:/CPRDData/Analysis_Murong/GAMLSS_imputed/Analysis 1/Model/Adult/", "m6_", infection, ".rds")    ##!! UPDATE ACCORDINGLY ##!!
    saveRDS(m6_infection, save_model_path)                                                                         ##!! UPDATE ACCORDINGLY ##!!
    logs <- c(logs, paste0("Model saved to: ", save_model_path))
    
    return(list(success = TRUE,    ##!! UPDATE ACCORDINGLY ##!!
                logs = logs))        
    
    
  },  error = function(e) {
    # Store the error message
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))  })
}



# (2) Set up parallel

# Infections for ADULTS
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti") #exclude "utio", "pyel", "pros", do them separately as they only includes one sex group

# exclude copd as it is very slow, do it seperately

num_cores <- 25

cl <- makeCluster(num_cores)

clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(gamlss))

clusterExport(cl, c("gamlss_adult", "infections"))

# (3) Run GAMLSS  

results <- parLapply(cl, infections, function(infection) gamlss_adult(infection, sample_size = 1))

stopCluster(cl)

beep("ping")





























