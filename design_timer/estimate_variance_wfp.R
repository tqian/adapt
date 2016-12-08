rm(list = ls())

setwd("~/Adaptive-Enrichment-Designs-with-Delayed-Outcomes/for_Josh/design_timer")


# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_delayWY_noL_wfp_delayWL0_delayLY*"
# file_prefix <- "gathered_estvar_AD47_delayWY_noL_wfp_delayWL0_delayLY" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_delayWY_noL_wfp.rds" # file to save
# ##### end of specification

##### specify the following #####
file_pattern <- "gathered_estvar_AD47_pW_wfp_pW*"
file_prefix <- "gathered_estvar_AD47_pW_wfp_pW" # for extracting param
file_suffix <- ".rda" # for extracting param
timer_filename <- "timers_AD47_pW_wfp.rds" # file to save
##### end of specification

# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_pW_grid_wfp*"
# file_prefix <- "gathered_estvar_AD47_pW_grid_wfp_case" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_pW_grid_wfp.rds" # file to save
# ##### end of specification

# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_pW_sameRsq_wfp*"
# file_prefix <- "gathered_estvar_AD47_pW_sameRsq_wfp_case" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_pW_sameRsq_wfp.rds" # file to save
# ##### end of specification

# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_pL_grid_wfp*"
# file_prefix <- "gathered_estvar_AD47_pL_grid_wfp_case" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_pL_grid_wfp.rds" # file to save
# ##### end of specification

# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_pL_sameRsq_wfp*"
# file_prefix <- "gathered_estvar_AD47_pL_sameRsq_wfp_case" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_pL_sameRsq_wfp.rds" # file to save
# ##### end of specification

# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_erate_dY4_wfp_erate*"
# file_prefix <- "gathered_estvar_AD47_erate_dY4_wfp_erate" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_erate_dY4_wfp.rds" # file to save
# ##### end of specification


# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_erate_dY1_wfp_erate*"
# file_prefix <- "gathered_estvar_AD47_erate_dY1_wfp_erate" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_erate_dY1_wfp.rds" # file to save
# ##### end of specification


# ##### specify the following #####
# file_pattern <- "gathered_estvar_AD47_pL_wfp_pL*"
# file_prefix <- "gathered_estvar_AD47_pL_wfp_pL" # for extracting param
# file_suffix <- ".rda" # for extracting param
# timer_filename <- "timers_AD47_pL_wfp.rds" # file to save
# ##### end of specification


# load("../results_simulation_paper_gathered/gathered_AD47_delayWY_noL_delayWL0_delayLY2.rda")
# vars <- diag(trials_H01H02$cov_unadj)
# vars <- matrix(vars, nrow = 3)
# vars_c <- vars[1, ]
# vars_1 <- vars[2, ]
# vars_2 <- vars[3, ]

vars_c <- c(0.036850023, 0.024491767, 0.019959962, 0.012448188, 0.007218001)
vars_1 <- c(0.07718200, 0.04946307, 0.04021796, 0.02495109, 0.01447511)
vars_2 <- c(0.07465222, 0.04823433, 0.03895192, 0.02434506, 0.01436538)
vars_c_inv <- vars_c^-1
vars_1_inv <- vars_1^-1
vars_2_inv <- vars_2^-1

K <- 5

setwd("../results_simulation_paper_estvar/")

files_all <- list.files()
files <- files_all[grepl(file_pattern, files_all)]
params <- sub(file_prefix, "", files)
params <- sub(file_suffix, "", params)
params <- as.numeric(params)

files <- files[order(params)]
params <- sort(params)

print(files)
print(params)

par(mfrow = c(2,2))

timer_design <- list()


find_single_timing <- function(desired_var_inv, var_inv_vector, timer_vector){
  right <- min(which(var_inv_vector > desired_var_inv))
  left <- right - 1
  
  v1 <- var_inv_vector[left]
  v2 <- var_inv_vector[right]
  v <- desired_var_inv
  t1 <- timer_vector[left]
  t2 <- timer_vector[right]
  
  t <- t2 - (t2-t1) * (v2-v) / (v2-v1)
  return(t)
}

for (est in c("ltmle", "unadj")) {
  
  timers <- data.frame(matrix(NA, nrow = length(files), ncol = K))
  rownames(timers) <- params
  colnames(timers) <- 1:K
  
  
  for (ifile in 1:length(files)){
    load(files[ifile])
    
    timer <- trial_info$timer
    vars <- diag(meanvar[[paste0("cov_", est)]])
    var_mat <- matrix(nrow = length(timer) * 2, ncol = 3)
    
    for (subpop in 0:2) {
      T_id <- seq(from = subpop*2 + 1, by = 6, length = length(timer))
      tT_id <- seq(from = subpop*2 + 2, by = 6, length = length(timer))
      v <- vars[T_id]
      tv <- vars[tT_id]
      vi <- v^-1
      tvi <- tv^-1
      
      var_mat[, subpop + 1] <- c(vi, tvi)
    }
    
    var_mat <- data.frame(var_mat)
    colnames(var_mat) <- c("comb_pop", "subpop1", "subpop2")
    var_mat$wait_for_pipeline <- rep(c("F", "T"), each = length(timer))
    var_mat$timer <- c(timer, timer)
    
    for (k in 1:K) {
      wfp <- ifelse(k == K, 'T', 'F')
      t_c <- find_single_timing(vars_c_inv[k], subset(var_mat, wait_for_pipeline == wfp)$comb_pop, timer)
      t_1 <- find_single_timing(vars_1_inv[k], subset(var_mat, wait_for_pipeline == wfp)$subpop1, timer)
      t_2 <- find_single_timing(vars_2_inv[k], subset(var_mat, wait_for_pipeline == wfp)$subpop2, timer)
      timers[ifile, k] <- max(t_c, t_1, t_2)
    }
  }
  
  timer_design[[est]] <- list(timers = timers)
}

print(timer_design)

setwd("../design_timer/")
saveRDS(timer_design, file = timer_filename)

