# Analysis 3: probability in antibiotic prescription
# This script run model specification 1 only (control for age sex ethnicity)

# Load necessary libraries
library(dplyr)
library(gamlss)
library(parallel)
library(doParallel)
library(beepr)
library(profvis)


################################################################################
# 1. Generate dataset and variable for analysis
################################################################################

# Data structure: one row per consultation, multiple rows per patient


# ---- Function to generate dataset and variables  ----

# Children 
generate_data_var_child <- function(infection) {
  
  tryCatch({
    message("Processing: ", infection)
    
    # ------------------            Data upload             -------------------#
    
    # Open disease dataset 
    input_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, ".rds") 
    infection_data <- readRDS(input_path)
    
    infection_data <- infection_data %>%   
      mutate(date = as.Date(date)) %>%
      arrange(patid, date)
    
    infection_data <- infection_data %>%
      dplyr::select(patid, abyes)
    
    # Open master dataset
    master_imputed <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds")
    
    master_select <- master_imputed %>%
      dplyr::select(patid, age, gender, imd, ethnic_2)
    
    
    # -----------   Generate child master disease visit dataset   -----------#
    
    # Merge datasets
    infection_sample_all <- infection_data %>%
      left_join(master_select, by = "patid")
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Child dataset ~~~~~~~~~~~~~~~~~~#
    
    infection_sample_child <- infection_sample_all %>%
      filter(age <18)
    
    # Generate B-spline for imd
    bs_imd <- bs(infection_sample_child$imd, df = 2, degree = 2)
    colnames(bs_imd) <- paste0("bs_imd_", seq_len(ncol(bs_imd)))
    attr_imd <- attributes(bs_imd)
    
    # Generate B-spline for age
    bs_age <- bs(infection_sample_child$age, df = 4)
    colnames(bs_age) <- paste0("bs_age_", seq_len(ncol(bs_age)))
    attr_age <- attributes(bs_age)
    
    # Combine with original data
    infection_sample_child <- cbind(infection_sample_child, bs_imd, bs_age)
    
    #Save to the new, safe location
    output_dir <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Data/Child/"
    
    output_path <- paste0(output_dir, infection, "_master_visit_child.rds") 
    output_path_imd <- paste0(output_dir, infection, "_visit_child_attr_imd.rds")
    output_path_age <- paste0(output_dir, infection, "_visit_child_attr_age.rds")
    
    saveRDS(infection_sample_child, output_path)
    saveRDS(attr_imd, output_path_imd)
    saveRDS(attr_age, output_path_age)
    
    return("Success")
  }, error = function(e) {
    return(paste("Error for", infection, ":", conditionMessage(e)))
  })
}




# Adult
generate_data_var_adult <- function(infection) {
  
  tryCatch({
    message("Processing: ", infection)
    
    # ------------------            Data upload             -------------------#
    
    # Open disease dataset 
    input_path <- paste0("D:/CPRDData/Analysis_Murong/Data/Disease_imputed/", infection, ".rds") 
    infection_data <- readRDS(input_path)
    
    infection_data <- infection_data %>%   
      mutate(date = as.Date(date)) %>%
      arrange(patid, date)
    
    infection_data <- infection_data %>%
      dplyr::select(patid, abyes)
    
    # Open master dataset
    master_imputed <- readRDS("D:/CPRDData/Analysis_Murong/Data/Multiple imputation/master_imputed.rds")
    
    master_select <- master_imputed %>%
      dplyr::select(patid, age, gender, imd, ethnic_2)
    
    
    # -----------   Generate adult master disease visit dataset   -----------#
    
    # Merge datasets
    infection_sample_all <- infection_data %>%
      left_join(master_select, by = "patid")
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Adult dataset ~~~~~~~~~~~~~~~~~~#
    
    infection_sample_adult <- infection_sample_all %>%
      filter(age >=18)
    
    # Generate B-spline for imd
    bs_imd <- bs(infection_sample_adult$imd, df = 2, degree = 2)
    colnames(bs_imd) <- paste0("bs_imd_", seq_len(ncol(bs_imd)))
    attr_imd <- attributes(bs_imd)
    
    # Generate B-spline for age
    bs_age <- bs(infection_sample_adult$age, df = 4)
    colnames(bs_age) <- paste0("bs_age_", seq_len(ncol(bs_age)))
    attr_age <- attributes(bs_age)
    
    # Combine with original data
    infection_sample_adult <- cbind(infection_sample_adult, bs_imd, bs_age)
    
    
    #Save to the new, safe location
    output_dir <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Data/Adult/"
    
    output_path <- paste0(output_dir, infection, "_master_visit_adult.rds") 
    output_path_imd <- paste0(output_dir, infection, "_visit_adult_attr_imd.rds")
    output_path_age <- paste0(output_dir, infection, "_visit_adult_attr_age.rds")
    
    saveRDS(infection_sample_adult, output_path)
    saveRDS(attr_imd, output_path_imd)
    saveRDS(attr_age, output_path_age)
    
    
    return("Success")
  }, error = function(e) {
    return(paste("Error for", infection, ":", conditionMessage(e)))
  })
}




# ---- Infections ----

#Infections for CHILDREN
infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", "utio", 
                "pyel", "sore", "sinu", "uurti", "bron", "ulrti")

# # Infections for ADULTS
infections <- c("asthma", "copd",  "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "gast", "utio", "pyel", "pros",
                "sore", "sinu", "uurti", "bron", "ulrti")

# --- Parallel Setup ---

num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Ensure the cluster is always stopped
on.exit(stopCluster(cl))


# --- Run in Parallel ---

# Note the addition of the .packages argument
results <- foreach(
  infection = infections,
  .packages = c("splines", "dplyr"),      
  .errorhandling = "pass"     
) %dopar% {
  generate_data_var_child(infection)   # !! UPDATE ACCORDINGLY !! _child/_adult #
}


names(results) <- infections
print("Processing summary:")
print(results)





################################################################################
# 2. Run model + calculate probability measure 
################################################################################

# Model: logistic regression + cluster SE on patient level 

# ----------    Children  -------------#

df_measure_child <- function(infection) {
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    # Paths
    data_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Data/Child/"
    measure_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Measure/Child/"
    
    
    
    # Load infection sample
    sample_file <- paste0(data_path, infection, "_master_visit_child.rds")   
    infection_sample <- readRDS(sample_file)
    assign("infection_sample", infection_sample, envir = .GlobalEnv)
    
    on.exit({
      rm(infection_sample, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # Load spline attributes
    attr_imd <- readRDS(paste0(data_path, infection, "_visit_child_attr_imd.rds"))
    attr_age <- readRDS(paste0(data_path, infection, "_visit_child_attr_age.rds"))
    logs$steps <- c(logs$steps, "Spline attributes loaded")
    
    
    # --------             Run model              ----------
    
    
    # Load model
    m1 <- glm(abyes ~ bs_imd_1 + bs_imd_2 +
                bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 +
                gender + ethnic_2,
              family = binomial(link = "logit"),
              data = infection_sample)
    
    logs$steps <- c(logs$steps, "Model completed")
    
    vcov_cluster <- sandwich::vcovCL(m1, cluster = ~ patid)
    coeftest(m1, vcov = vcov_cluster)
    
    logs$steps <- c(logs$steps, "Cluster SE calculated")
    
    
    # -------- Create new data set for prediction  ----------
    
    # Spline basis for IMD 
    imd_seq <- seq(min(infection_sample$imd), max(infection_sample$imd), length.out = 100)
    bs_imd_seq <- bs(imd_seq, knots = attr_imd$knots, degree = attr_imd$degree, Boundary.knots = attr_imd$Boundary.knots)
    colnames(bs_imd_seq) <- paste0("bs_imd_", seq_len(ncol(bs_imd_seq)))
    logs$steps <- c(logs$steps, "IMD spline basis created")
    
    
    # Spline basis for constant age
    age_ref <- median(infection_sample$age)
    bs_age_seq <- bs(age_ref, knots = attr_age$knots, degree = attr_age$degree, Boundary.knots = attr_age$Boundary.knots)
    colnames(bs_age_seq) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created")
    
    n_row <- nrow(bs_imd_seq)
    bs_age_repeated <- matrix(rep(bs_age_seq, each = n_row), ncol = ncol(bs_age_seq), byrow = FALSE)
    colnames(bs_age_repeated) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created and repeated")
    
    
    # Reference values
    ref_vals <- data.frame(
      age = age_ref,
      gender = factor(levels(infection_sample$gender)[1], levels = levels(infection_sample$gender)),
      ethnic_2 = factor(levels(infection_sample$ethnic_2)[1], levels = levels(infection_sample$ethnic_2))
    )
    logs$steps <- c(logs$steps, "Reference values created")
    
    # New data
    newdata <- cbind(ref_vals[rep(1, n_row), ], bs_imd_seq, bs_age_repeated)
    logs$steps <- c(logs$steps, "New data for prediction created")
    
    
    # -------- Simulation ----------
    
    coef_all <- coef(m1)
    logs$steps <- c(logs$steps, "coef_all created")
    
    vcov_all <- vcov_cluster
    logs$steps <- c(logs$steps, "vcov_all created")
    
    
    sim_coefs <- MASS::mvrnorm(1000, coef_all, vcov_all)
    logs$steps <- c(logs$steps, "Simulated 1000 sets of coefficients")
    
    Xp_mu <- model.matrix(formula(m1, "mu")[-2], data = newdata)
    
    idx_mu <- 1:ncol(Xp_mu)
    mu_sim <- Xp_mu %*% t(sim_coefs[, idx_mu])
    logs$steps <- c(logs$steps, "Predictions simulated")
    
    mu_resp <- 1 / (1 + exp(-(mu_sim)))
    logs$steps <- c(logs$steps, "Response variables calculated")
    
    
    # -------- summaries mean and variance ----------
    
    summarise_sim <- function(mat) {
      apply(mat, 1, function(x) c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975)))
    }
    df_mean <- as.data.frame(t(summarise_sim(mu_resp)))
    df_mean$imd <- imd_seq
    
    names(df_mean)[2:3] <- c("lower","upper")
    
    logs$steps <- c(logs$steps, "Mean and variance calculated")
    
    # -------- save mean and variance ----------
    
    saveRDS(df_mean, paste0(measure_path,"df_mean_", infection, "_m1.rds"))
    
    logs$steps <- c(logs$steps, "Mean df saved")
    
    
    logs$end_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    
    return(list(success = TRUE, logs = logs))
    
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))
  })
  
}

# Run parallel ###

infections <- c("asthma", "lrti", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", 
                "sore", "sinu", "uurti", "bron", "ulrti")

num_cores <- 20

cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(dplyr)
  library(splines)
  library(sandwich)
  library(lmtest)
})


clusterExport(cl, c("df_measure_child", "infections"))

# (3) Parallel  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- parLapply(cl, infections, df_measure_child) 

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")

print(results)


# UTI


df_measure_child <- function(infection) {
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    # Paths
    data_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Data/Child/"
    measure_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Measure/Child/"
    
    # Load infection sample
    sample_file <- paste0(data_path, infection, "_master_visit_child.rds")   
    infection_sample <- readRDS(sample_file)
    assign("infection_sample", infection_sample, envir = .GlobalEnv)
    
    on.exit({
      rm(infection_sample, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # Load spline attributes
    attr_imd <- readRDS(paste0(data_path, infection, "_visit_child_attr_imd.rds"))
    attr_age <- readRDS(paste0(data_path, infection, "_visit_child_attr_age.rds"))
    logs$steps <- c(logs$steps, "Spline attributes loaded")
    
    
    # --------             Run model              ----------
    
    
    # Load model
    m1 <- glm(abyes ~ bs_imd_1 + bs_imd_2 +
                bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + ethnic_2,
              family = binomial(link = "logit"),
              data = infection_sample)
    
    logs$steps <- c(logs$steps, "Model completed")
    
    vcov_cluster <- sandwich::vcovCL(m1, cluster = ~ patid)
    coeftest(m1, vcov = vcov_cluster)
    
    logs$steps <- c(logs$steps, "Cluster SE calculated")
    
    
    # -------- Create new data set for prediction  ----------
    
    # Spline basis for IMD 
    imd_seq <- seq(min(infection_sample$imd), max(infection_sample$imd), length.out = 100)
    bs_imd_seq <- bs(imd_seq, knots = attr_imd$knots, degree = attr_imd$degree, Boundary.knots = attr_imd$Boundary.knots)
    colnames(bs_imd_seq) <- paste0("bs_imd_", seq_len(ncol(bs_imd_seq)))
    logs$steps <- c(logs$steps, "IMD spline basis created")
    
    
    # Spline basis for constant age
    age_ref <- median(infection_sample$age)
    bs_age_seq <- bs(age_ref, knots = attr_age$knots, degree = attr_age$degree, Boundary.knots = attr_age$Boundary.knots)
    colnames(bs_age_seq) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created")
    
    n_row <- nrow(bs_imd_seq)
    bs_age_repeated <- matrix(rep(bs_age_seq, each = n_row), ncol = ncol(bs_age_seq), byrow = FALSE)
    colnames(bs_age_repeated) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created and repeated")
    
    
    # Reference values
    ref_vals <- data.frame(
      age = age_ref,
      ethnic_2 = factor(levels(infection_sample$ethnic_2)[1], levels = levels(infection_sample$ethnic_2))
    )
    logs$steps <- c(logs$steps, "Reference values created")
    
    # New data
    newdata <- cbind(ref_vals[rep(1, n_row), ], bs_imd_seq, bs_age_repeated)
    logs$steps <- c(logs$steps, "New data for prediction created")
    
    
    # -------- Simulation ----------
    
    coef_all <- coef(m1)
    logs$steps <- c(logs$steps, "coef_all created")
    
    vcov_all <- vcov_cluster
    logs$steps <- c(logs$steps, "vcov_all created")
    
    
    sim_coefs <- MASS::mvrnorm(1000, coef_all, vcov_all)
    logs$steps <- c(logs$steps, "Simulated 1000 sets of coefficients")
    
    Xp_mu <- model.matrix(formula(m1, "mu")[-2], data = newdata)
    
    idx_mu <- 1:ncol(Xp_mu)
    mu_sim <- Xp_mu %*% t(sim_coefs[, idx_mu])
    logs$steps <- c(logs$steps, "Predictions simulated")
    
    mu_resp <- 1 / (1 + exp(-(mu_sim)))
    logs$steps <- c(logs$steps, "Response variables calculated")
    
    
    # -------- summaries mean and variance ----------
    
    summarise_sim <- function(mat) {
      apply(mat, 1, function(x) c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975)))
    }
    df_mean <- as.data.frame(t(summarise_sim(mu_resp)))
    df_mean$imd <- imd_seq
    
    names(df_mean)[2:3] <- c("lower","upper")
    
    logs$steps <- c(logs$steps, "Mean and variance calculated")
    
    # -------- save mean and variance ----------
    
    saveRDS(df_mean, paste0(measure_path,"df_mean_", infection, "_m1.rds"))
    
    logs$steps <- c(logs$steps, "Mean df saved")
    
    
    logs$end_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    
    return(list(success = TRUE, logs = logs))
    
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))
  })
  
}

# Run parallel ###

infections <- c("utio")

num_cores <- 5

cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(dplyr)
  library(splines)
  library(sandwich)
  library(lmtest)
})


clusterExport(cl, c("df_measure_child", "infections"))

# (3) Run GAMLSS  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- parLapply(cl, infections, df_measure_child) 

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")

print(results)




# ----------    Adults  -------------#


df_measure_adult <- function(infection) {
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    # Paths
    data_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Data/Adult/"
    measure_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Measure/Adult/"
    
    # Load infection sample
    sample_file <- paste0(data_path, infection, "_master_visit_adult.rds")   
    infection_sample <- readRDS(sample_file)
    assign("infection_sample", infection_sample, envir = .GlobalEnv)
    
    on.exit({
      rm(infection_sample, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # Load spline attributes
    attr_imd <- readRDS(paste0(data_path, infection, "_visit_adult_attr_imd.rds"))
    attr_age <- readRDS(paste0(data_path, infection, "_visit_adult_attr_age.rds"))
    logs$steps <- c(logs$steps, "Spline attributes loaded")
    
    
    # --------             Run model              ----------
    
    
    # Load model
    m1 <- glm(abyes ~ bs_imd_1 + bs_imd_2 +
                bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 +
                gender + ethnic_2,
              family = binomial(link = "logit"),
              data = infection_sample)
    
    logs$steps <- c(logs$steps, "Model completed")
    
    vcov_cluster <- sandwich::vcovCL(m1, cluster = ~ patid)
    coeftest(m1, vcov = vcov_cluster)
    
    logs$steps <- c(logs$steps, "Cluster SE calculated")
    
    
    # -------- Create new data set for prediction  ----------
    
    # Spline basis for IMD 
    imd_seq <- seq(min(infection_sample$imd), max(infection_sample$imd), length.out = 100)
    bs_imd_seq <- bs(imd_seq, knots = attr_imd$knots, degree = attr_imd$degree, Boundary.knots = attr_imd$Boundary.knots)
    colnames(bs_imd_seq) <- paste0("bs_imd_", seq_len(ncol(bs_imd_seq)))
    logs$steps <- c(logs$steps, "IMD spline basis created")
    
    
    # Spline basis for constant age
    age_ref <- median(infection_sample$age)
    bs_age_seq <- bs(age_ref, knots = attr_age$knots, degree = attr_age$degree, Boundary.knots = attr_age$Boundary.knots)
    colnames(bs_age_seq) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created")
    
    n_row <- nrow(bs_imd_seq)
    bs_age_repeated <- matrix(rep(bs_age_seq, each = n_row), ncol = ncol(bs_age_seq), byrow = FALSE)
    colnames(bs_age_repeated) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created and repeated")
    
    
    # Reference values
    ref_vals <- data.frame(
      age = age_ref,
      gender = factor(levels(infection_sample$gender)[1], levels = levels(infection_sample$gender)),
      ethnic_2 = factor(levels(infection_sample$ethnic_2)[1], levels = levels(infection_sample$ethnic_2))
    )
    logs$steps <- c(logs$steps, "Reference values created")
    
    # New data
    newdata <- cbind(ref_vals[rep(1, n_row), ], bs_imd_seq, bs_age_repeated)
    logs$steps <- c(logs$steps, "New data for prediction created")
    
    
    # -------- Simulation ----------
    
    coef_all <- coef(m1)
    logs$steps <- c(logs$steps, "coef_all created")
    
    vcov_all <- vcov_cluster
    logs$steps <- c(logs$steps, "vcov_all created")
    
    
    sim_coefs <- MASS::mvrnorm(1000, coef_all, vcov_all)
    logs$steps <- c(logs$steps, "Simulated 1000 sets of coefficients")
    
    Xp_mu <- model.matrix(formula(m1, "mu")[-2], data = newdata)
    
    idx_mu <- 1:ncol(Xp_mu)
    mu_sim <- Xp_mu %*% t(sim_coefs[, idx_mu])
    logs$steps <- c(logs$steps, "Predictions simulated")
    
    mu_resp <- 1 / (1 + exp(-(mu_sim)))
    logs$steps <- c(logs$steps, "Response variables calculated")
    
    
    # -------- summaries mean and variance ----------
    
    summarise_sim <- function(mat) {
      apply(mat, 1, function(x) c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975)))
    }
    df_mean <- as.data.frame(t(summarise_sim(mu_resp)))
    df_mean$imd <- imd_seq
    
    names(df_mean)[2:3] <- c("lower","upper")
    
    logs$steps <- c(logs$steps, "Mean and variance calculated")
    
    # -------- save mean and variance ----------
    
    saveRDS(df_mean, paste0(measure_path,"df_mean_", infection, "_m1.rds"))
    
    logs$steps <- c(logs$steps, "Mean df saved")
    
    
    logs$end_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    
    return(list(success = TRUE, logs = logs))
    
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))
  })
  
}

# Run parallel ###

infections <-   c("asthma", "lrti", "copd", "cough", "urti","pneu", "cell", "impe", "exter", "media", "gast", "utio", "pyel", 
                  "sore", "sinu", "uurti", "bron", "ulrti") #except for "utio", "pyel", "pros"

num_cores <- 20

cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(dplyr)
  library(splines)
  library(sandwich)
  library(lmtest)
})


clusterExport(cl, c("df_measure_adult", "infections"))

# (3) Run GAMLSS  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- parLapply(cl, infections, df_measure_adult) 

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")

print(results)



# uti, pros, pyel 


df_measure_adult <- function(infection) {
  logs <- list()
  logs$start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logs$infection <- infection
  logs$steps <- c()
  
  logs$steps <- c(logs$steps, paste("Running infection:", infection))
  
  tryCatch({
    # Paths
    data_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Data/Adult/"
    measure_path <- "D:/CPRDData/Analysis_Murong/GAMLSS_imputed_observed/Analysis 3/Measure/Adult/"
    
    # Load infection sample
    sample_file <- paste0(data_path, infection, "_master_visit_adult.rds")   
    infection_sample <- readRDS(sample_file)
    assign("infection_sample", infection_sample, envir = .GlobalEnv)
    
    on.exit({
      rm(infection_sample, envir = .GlobalEnv)
      gc()
    }, add = TRUE)
    
    logs$steps <- c(logs$steps, "Uploaded data")
    
    # Load spline attributes
    attr_imd <- readRDS(paste0(data_path, infection, "_visit_adult_attr_imd.rds"))
    attr_age <- readRDS(paste0(data_path, infection, "_visit_adult_attr_age.rds"))
    logs$steps <- c(logs$steps, "Spline attributes loaded")
    
    
    # --------             Run model              ----------
    
    
    # Load model
    m1 <- glm(abyes ~ bs_imd_1 + bs_imd_2 +
                bs_age_1 + bs_age_2 + bs_age_3 + bs_age_4 + ethnic_2,
              family = binomial(link = "logit"),
              data = infection_sample)
    
    logs$steps <- c(logs$steps, "Model completed")
    
    vcov_cluster <- sandwich::vcovCL(m1, cluster = ~ patid)
    coeftest(m1, vcov = vcov_cluster)
    
    logs$steps <- c(logs$steps, "Cluster SE calculated")
    
    
    # -------- Create new data set for prediction  ----------
    
    # Spline basis for IMD 
    imd_seq <- seq(min(infection_sample$imd), max(infection_sample$imd), length.out = 100)
    bs_imd_seq <- bs(imd_seq, knots = attr_imd$knots, degree = attr_imd$degree, Boundary.knots = attr_imd$Boundary.knots)
    colnames(bs_imd_seq) <- paste0("bs_imd_", seq_len(ncol(bs_imd_seq)))
    logs$steps <- c(logs$steps, "IMD spline basis created")
    
    
    # Spline basis for constant age
    age_ref <- median(infection_sample$age)
    bs_age_seq <- bs(age_ref, knots = attr_age$knots, degree = attr_age$degree, Boundary.knots = attr_age$Boundary.knots)
    colnames(bs_age_seq) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created")
    
    n_row <- nrow(bs_imd_seq)
    bs_age_repeated <- matrix(rep(bs_age_seq, each = n_row), ncol = ncol(bs_age_seq), byrow = FALSE)
    colnames(bs_age_repeated) <- paste0("bs_age_", seq_len(ncol(bs_age_seq)))
    logs$steps <- c(logs$steps, "Age spline basis created and repeated")
    
    
    # Reference values
    ref_vals <- data.frame(
      age = age_ref,
      ethnic_2 = factor(levels(infection_sample$ethnic_2)[1], levels = levels(infection_sample$ethnic_2))
    )
    logs$steps <- c(logs$steps, "Reference values created")
    
    # New data
    newdata <- cbind(ref_vals[rep(1, n_row), ], bs_imd_seq, bs_age_repeated)
    logs$steps <- c(logs$steps, "New data for prediction created")
    
    
    # -------- Simulation ----------
    
    coef_all <- coef(m1)
    logs$steps <- c(logs$steps, "coef_all created")
    
    vcov_all <- vcov_cluster
    logs$steps <- c(logs$steps, "vcov_all created")
    
    
    sim_coefs <- MASS::mvrnorm(1000, coef_all, vcov_all)
    logs$steps <- c(logs$steps, "Simulated 1000 sets of coefficients")
    
    Xp_mu <- model.matrix(formula(m1, "mu")[-2], data = newdata)
    
    idx_mu <- 1:ncol(Xp_mu)
    mu_sim <- Xp_mu %*% t(sim_coefs[, idx_mu])
    logs$steps <- c(logs$steps, "Predictions simulated")
    
    mu_resp <- 1 / (1 + exp(-(mu_sim)))
    logs$steps <- c(logs$steps, "Response variables calculated")
    
    
    # -------- summaries mean and variance ----------
    
    summarise_sim <- function(mat) {
      apply(mat, 1, function(x) c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975)))
    }
    df_mean <- as.data.frame(t(summarise_sim(mu_resp)))
    df_mean$imd <- imd_seq
    
    names(df_mean)[2:3] <- c("lower","upper")
    
    logs$steps <- c(logs$steps, "Mean and variance calculated")
    
    # -------- save mean and variance ----------
    
    saveRDS(df_mean, paste0(measure_path,"df_mean_", infection, "_m1.rds"))
    
    logs$steps <- c(logs$steps, "Mean df saved")
    
    
    logs$end_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    
    return(list(success = TRUE, logs = logs))
    
    
  }, error = function(e){
    logs$error <- conditionMessage(e)
    logs$trace <- traceback()
    return(list(success = FALSE, logs = logs))
  })
  
}

# Run parallel ###

infections <-   c("utio", "pyel", "pros")

num_cores <- 5

cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(dplyr)
  library(splines)
  library(sandwich)
  library(lmtest)
})


clusterExport(cl, c("df_measure_adult", "infections"))

# (3) Run GAMLSS  

cat("Full process/parallel execution started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- parLapply(cl, infections, df_measure_adult) 

stopCluster(cl)

cat("Parallel execution ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

beep("ping")

print(results)
