# Gather all the parallel computing results from simtrials, then
# compute the covariance matrix and correlation matrix, for ltmle estimators and unadj estimators, respectively.

# Output: a list of: 
# 1) 4 matrices of all the estimators: ltmle & unadj, original & standardized (Z-statistic);
# 2) 4 covariance matrices of the above 4 matrices of estimators;
# 3) 4 mean vectors of the above 4 matices of estimators


# for simu_eval estimate variance on the fly within each trial:

gather_stack_result <- function(filenames) {
  stack_table_byeffset <- lapply(filenames, function(names) {
    stacktable <- readRDS(names[1])
    
    if (length(names >= 2)) {
      for (i in 2:length(names)) {
        stacktable <- rbind(stacktable, readRDS(names[i]))
      }
    }
    return(stacktable)
  })
  return(stack_table_byeffset)
}





gen_names <- function(generic_filename, idx){
  # generic_filename: the part of the filename that is shared across all the result files.
  # idx: the indices of the parallel jobs (for example, 1:50)
  
  names <- paste0(generic_filename, "_", idx, ".rda")
  return(names)
}


gather_results_timer <- function(filenames, wait_for_pipeline = FALSE){
  
  n <- length(filenames)
  
  load(filenames[1])
  K <- ncol(test$result$subpop_c$mat_ltmle)
  
  num_subpop <- 2
  
  ncol <- ifelse(wait_for_pipeline, 2*K*(num_subpop+1), K*(num_subpop+1))
  
  # initialize overall matrix to collect all parallel simulated esitmators
  mat_ltmle <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj <- matrix(NA, nrow = 1, ncol = ncol)
  mat_ltmle_a1 <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj_a1 <- matrix(NA, nrow = 1, ncol = ncol)
  mat_ltmle_a0 <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj_a0 <- matrix(NA, nrow = 1, ncol = ncol)
  
  mat_ltmle_ICvar <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj_samplevar <- matrix(NA, nrow = 1, ncol = ncol)
  
  ss_enrolled_c <- ss_L_observed_c <- ss_Y_observed_c <- matrix(NA, nrow = 1, ncol = K)
  ss_enrolled_1 <- ss_L_observed_1 <- ss_Y_observed_1 <- matrix(NA, nrow = 1, ncol = K)
  ss_enrolled_2 <- ss_L_observed_2 <- ss_Y_observed_2 <- matrix(NA, nrow = 1, ncol = K)
  
  # rbind each resulting matrix to the overall matrix  
  for (i in 1:n){
    load(filenames[i])
    tmp_ltmle <- cbind(test$result$subpop_c$mat_ltmle, test$result$subpop_1$mat_ltmle, test$result$subpop_2$mat_ltmle)
    mat_ltmle <- rbind(mat_ltmle, tmp_ltmle)
    tmp_unadj <- cbind(test$result$subpop_c$mat_unadj, test$result$subpop_1$mat_unadj, test$result$subpop_2$mat_unadj)
    mat_unadj <- rbind(mat_unadj, tmp_unadj)
    
    tmp_ltmle_a1 <- cbind(test$result_a1$subpop_c$mat_ltmle, test$result_a1$subpop_1$mat_ltmle, test$result_a1$subpop_2$mat_ltmle)
    mat_ltmle_a1 <- rbind(mat_ltmle_a1, tmp_ltmle_a1)
    tmp_unadj_a1 <- cbind(test$result_a1$subpop_c$mat_unadj, test$result_a1$subpop_1$mat_unadj, test$result_a1$subpop_2$mat_unadj)
    mat_unadj_a1 <- rbind(mat_unadj_a1, tmp_unadj_a1)
    
    tmp_ltmle_a0 <- cbind(test$result_a0$subpop_c$mat_ltmle, test$result_a0$subpop_1$mat_ltmle, test$result_a0$subpop_2$mat_ltmle)
    mat_ltmle_a0 <- rbind(mat_ltmle_a0, tmp_ltmle_a0)
    tmp_unadj_a0 <- cbind(test$result_a0$subpop_c$mat_unadj, test$result_a0$subpop_1$mat_unadj, test$result_a0$subpop_2$mat_unadj)
    mat_unadj_a0 <- rbind(mat_unadj_a0, tmp_unadj_a0)
    
    if ("mat_ltmle_ICvar" %in% names(test$result$subpop_c)) {
      tmp_ltmle_ICvar <- cbind(test$result$subpop_c$mat_ltmle_ICvar, test$result$subpop_1$mat_ltmle_ICvar, test$result$subpop_2$mat_ltmle_ICvar)
      mat_ltmle_ICvar <- rbind(mat_ltmle_ICvar, tmp_ltmle_ICvar)
    }
    
    if ("mat_unadj_samplevar" %in% names(test$result$subpop_c)) {
      tmp_unadj_samplevar <- cbind(test$result$subpop_c$mat_unadj_samplevar, test$result$subpop_1$mat_unadj_samplevar, test$result$subpop_2$mat_unadj_samplevar)
      mat_unadj_samplevar <- rbind(mat_unadj_samplevar, tmp_unadj_samplevar)
    }
    
    ss_enrolled_c <- rbind(ss_enrolled_c, test$ss_enrolled_c)
    ss_L_observed_c <- rbind(ss_L_observed_c, test$ss_L_observed_c)
    ss_Y_observed_c <- rbind(ss_Y_observed_c, test$ss_Y_observed_c)
    
    ss_enrolled_1 <- rbind(ss_enrolled_1, test$ss_enrolled_1)
    ss_L_observed_1 <- rbind(ss_L_observed_1, test$ss_L_observed_1)
    ss_Y_observed_1 <- rbind(ss_Y_observed_1, test$ss_Y_observed_1)
    
    ss_enrolled_2 <- rbind(ss_enrolled_2, test$ss_enrolled_2)
    ss_L_observed_2 <- rbind(ss_L_observed_2, test$ss_L_observed_2)
    ss_Y_observed_2 <- rbind(ss_Y_observed_2, test$ss_Y_observed_2)
    ###### Remove the ss first rows
  }
  
  # delete the first line (initializer)
  mat_ltmle <- mat_ltmle[-1, ]
  mat_unadj <- mat_unadj[-1, ]
  mat_ltmle_a1 <- mat_ltmle_a1[-1, ]
  mat_unadj_a1 <- mat_unadj_a1[-1, ]
  mat_ltmle_a0 <- mat_ltmle_a0[-1, ]
  mat_unadj_a0 <- mat_unadj_a0[-1, ]
  mat_ltmle_ICvar <- mat_ltmle_ICvar[-1, ]
  mat_unadj_samplevar <- mat_unadj_samplevar[-1, ]
  
  ss_enrolled_c <- ss_enrolled_c[-1, ]
  ss_L_observed_c <- ss_L_observed_c[-1, ]
  ss_Y_observed_c <- ss_Y_observed_c[-1, ]
  
  ss_enrolled_1 <- ss_enrolled_1[-1, ]
  ss_L_observed_1 <- ss_L_observed_1[-1, ]
  ss_Y_observed_1 <- ss_Y_observed_1[-1, ]
  
  ss_enrolled_2 <- ss_enrolled_2[-1, ]
  ss_L_observed_2 <- ss_L_observed_2[-1, ]
  ss_Y_observed_2 <- ss_Y_observed_2[-1, ]
  
  # swap columns
  if (wait_for_pipeline) {
    curr_id <- as.vector(t(outer(c("Tc_", "tTc_", "T1_", "tT1_", "T2_", "tT2_"), 1:K, FUN = paste0)))
    swap_id <- as.vector(outer(c("Tc_", "tTc_", "T1_", "tT1_", "T2_", "tT2_"), 1:K, FUN = paste0))
  } else {
    curr_id <- as.vector(t(outer(c("Tc_", "T1_", "T2_"), 1:K, FUN = paste0)))
    swap_id <- as.vector(outer(c("Tc_", "T1_", "T2_"), 1:K, FUN = paste0))
  }
  
  colnames(mat_ltmle) <- colnames(mat_unadj) <-
    colnames(mat_ltmle_a1) <- colnames(mat_unadj_a1) <-
    colnames(mat_ltmle_a0) <- colnames(mat_unadj_a0) <-
    colnames(mat_ltmle_ICvar) <- colnames(mat_unadj_samplevar) <- curr_id
  
  mat_ltmle <- mat_ltmle[, swap_id]
  mat_unadj <- mat_unadj[, swap_id]
  mat_ltmle_a1 <- mat_ltmle_a1[, swap_id]
  mat_unadj_a1 <- mat_unadj_a1[, swap_id]
  mat_ltmle_a0 <- mat_ltmle_a0[, swap_id]
  mat_unadj_a0 <- mat_unadj_a0[, swap_id]
  mat_ltmle_ICvar <- mat_ltmle_ICvar[, swap_id]
  mat_unadj_samplevar <- mat_unadj_samplevar[, swap_id]
  
  # Compute mean and covariance of ltmle and unadj
  mean_ltmle <- apply(mat_ltmle, 2, mean)
  cov_ltmle <- cov(mat_ltmle)
  mean_unadj <- apply(mat_unadj, 2, mean)
  cov_unadj <- cov(mat_unadj)
  
  mean_ltmle_a1 <- apply(mat_ltmle_a1, 2, mean)
  cov_ltmle_a1 <- cov(mat_ltmle_a1)
  mean_unadj_a1 <- apply(mat_unadj_a1, 2, mean)
  cov_unadj_a1 <- cov(mat_unadj_a1)
  
  mean_ltmle_a0 <- apply(mat_ltmle_a0, 2, mean)
  cov_ltmle_a0 <- cov(mat_ltmle_a0)
  mean_unadj_a0 <- apply(mat_unadj_a0, 2, mean)
  cov_unadj_a0 <- cov(mat_unadj_a0)
  
  # Standardize ltmle and unadj (divide by standard error, make it Z-statisitc)
  mat_ltmle_std <- mat_ltmle %*% diag(diag(cov_ltmle) ^ (-1/2))
  mat_unadj_std <- mat_unadj %*% diag(diag(cov_unadj) ^ (-1/2))
  
  mat_ltmle_a1_std <- mat_ltmle_a1 %*% diag(diag(cov_ltmle_a1) ^ (-1/2))
  mat_unadj_a1_std <- mat_unadj_a1 %*% diag(diag(cov_unadj_a1) ^ (-1/2))
  
  mat_ltmle_a0_std <- mat_ltmle_a0 %*% diag(diag(cov_ltmle_a0) ^ (-1/2))
  mat_unadj_a0_std <- mat_unadj_a0 %*% diag(diag(cov_unadj_a0) ^ (-1/2))
  
  # Compute mean and covariance of the standardized ones
  mean_ltmle_std <- apply(mat_ltmle_std, 2, mean)
  cov_ltmle_std <- cov(mat_ltmle_std)
  mean_unadj_std <- apply(mat_unadj_std, 2, mean)
  cov_unadj_std <- cov(mat_unadj_std)
  
  mean_ltmle_a1_std <- apply(mat_ltmle_a1_std, 2, mean)
  cov_ltmle_a1_std <- cov(mat_ltmle_a1_std)
  mean_unadj_a1_std <- apply(mat_unadj_a1_std, 2, mean)
  cov_unadj_a1_std <- cov(mat_unadj_a1_std)
  
  mean_ltmle_a0_std <- apply(mat_ltmle_a0_std, 2, mean)
  cov_ltmle_a0_std <- cov(mat_ltmle_a0_std)
  mean_unadj_a0_std <- apply(mat_unadj_a0_std, 2, mean)
  cov_unadj_a0_std <- cov(mat_unadj_a0_std)
  
  # gather the cov/corr results
  result <- list(mat_ltmle = mat_ltmle, mat_unadj = mat_unadj,
                 mean_ltmle = mean_ltmle, mean_unadj = mean_unadj,
                 cov_ltmle = cov_ltmle, cov_unadj = cov_unadj,
                 mat_ltmle_std = mat_ltmle_std, mat_unadj_std = mat_unadj_std,
                 mean_ltmle_std = mean_ltmle_std, mean_unadj_std = mean_unadj_std,
                 cov_ltmle_std = cov_ltmle_std, cov_unadj_std = cov_unadj_std,
                 mat_ltmle_ICvar = mat_ltmle_ICvar,
                 mat_unadj_samplevar = mat_unadj_samplevar,
                 
                 mat_ltmle_a1 = mat_ltmle_a1, mat_unadj_a1 = mat_unadj_a1,
                 mean_ltmle_a1 = mean_ltmle_a1, mean_unadj_a1 = mean_unadj_a1,
                 cov_ltmle_a1 = cov_ltmle_a1, cov_unadj_a1 = cov_unadj_a1,
                 mat_ltmle_a1_std = mat_ltmle_a1_std, mat_unadj_a1_std = mat_unadj_a1_std,
                 mean_ltmle_a1_std = mean_ltmle_a1_std, mean_unadj_a1_std = mean_unadj_a1_std,
                 cov_ltmle_a1_std = cov_ltmle_a1_std, cov_unadj_a1_std = cov_unadj_a1_std,
                 
                 mat_ltmle_a0 = mat_ltmle_a0, mat_unadj_a0 = mat_unadj_a0,
                 mean_ltmle_a0 = mean_ltmle_a0, mean_unadj_a0 = mean_unadj_a0,
                 cov_ltmle_a0 = cov_ltmle_a0, cov_unadj_a0 = cov_unadj_a0,
                 mat_ltmle_a0_std = mat_ltmle_a0_std, mat_unadj_a0_std = mat_unadj_a0_std,
                 mean_ltmle_a0_std = mean_ltmle_a0_std, mean_unadj_a0_std = mean_unadj_a0_std,
                 cov_ltmle_a0_std = cov_ltmle_a0_std, cov_unadj_a0_std = cov_unadj_a0_std,
                 
                 ss_enrolled_c = ss_enrolled_c, ss_L_observed_c = ss_L_observed_c, ss_Y_observed_c = ss_Y_observed_c,
                 ss_enrolled_1 = ss_enrolled_1, ss_L_observed_1 = ss_L_observed_1, ss_Y_observed_1 = ss_Y_observed_1,
                 ss_enrolled_2 = ss_enrolled_2, ss_L_observed_2 = ss_L_observed_2, ss_Y_observed_2 = ss_Y_observed_2)
  return(result)
}


gather_results <- function(filenames, wait_for_pipeline = FALSE){
  
  n <- length(filenames)
  
  load(filenames[1])
  
  K <- ifelse(wait_for_pipeline, ncol(test$result$subpop_c$mat_ltmle) / 2, 
              ncol(test$result$subpop_c$mat_ltmle))

  num_subpop <- 2
  
  ncol <- ncol(test$result$subpop_c$mat_ltmle) * (num_subpop + 1)
  
  # initialize overall matrix to collect all parallel simulated esitmators
  mat_ltmle <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj <- matrix(NA, nrow = 1, ncol = ncol)
  mat_ltmle_a1 <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj_a1 <- matrix(NA, nrow = 1, ncol = ncol)
  mat_ltmle_a0 <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj_a0 <- matrix(NA, nrow = 1, ncol = ncol)
  
  mat_ltmle_ICvar <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj_samplevar <- matrix(NA, nrow = 1, ncol = ncol)
  
  # rbind each resulting matrix to the overall matrix  
  for (i in 1:n){
    load(filenames[i])
    tmp_ltmle <- cbind(test$result$subpop_c$mat_ltmle, test$result$subpop_1$mat_ltmle, test$result$subpop_2$mat_ltmle)
    mat_ltmle <- rbind(mat_ltmle, tmp_ltmle)
    tmp_unadj <- cbind(test$result$subpop_c$mat_unadj, test$result$subpop_1$mat_unadj, test$result$subpop_2$mat_unadj)
    mat_unadj <- rbind(mat_unadj, tmp_unadj)
    
    tmp_ltmle_a1 <- cbind(test$result_a1$subpop_c$mat_ltmle, test$result_a1$subpop_1$mat_ltmle, test$result_a1$subpop_2$mat_ltmle)
    mat_ltmle_a1 <- rbind(mat_ltmle_a1, tmp_ltmle_a1)
    tmp_unadj_a1 <- cbind(test$result_a1$subpop_c$mat_unadj, test$result_a1$subpop_1$mat_unadj, test$result_a1$subpop_2$mat_unadj)
    mat_unadj_a1 <- rbind(mat_unadj_a1, tmp_unadj_a1)
    
    tmp_ltmle_a0 <- cbind(test$result_a0$subpop_c$mat_ltmle, test$result_a0$subpop_1$mat_ltmle, test$result_a0$subpop_2$mat_ltmle)
    mat_ltmle_a0 <- rbind(mat_ltmle_a0, tmp_ltmle_a0)
    tmp_unadj_a0 <- cbind(test$result_a0$subpop_c$mat_unadj, test$result_a0$subpop_1$mat_unadj, test$result_a0$subpop_2$mat_unadj)
    mat_unadj_a0 <- rbind(mat_unadj_a0, tmp_unadj_a0)
    
    if ("mat_ltmle_ICvar" %in% names(test$result$subpop_c)) {
      tmp_ltmle_ICvar <- cbind(test$result$subpop_c$mat_ltmle_ICvar, test$result$subpop_1$mat_ltmle_ICvar, test$result$subpop_2$mat_ltmle_ICvar)
      mat_ltmle_ICvar <- rbind(mat_ltmle_ICvar, tmp_ltmle_ICvar)
    }
    
    if ("mat_unadj_samplevar" %in% names(test$result$subpop_c)) {
      tmp_unadj_samplevar <- cbind(test$result$subpop_c$mat_unadj_samplevar, test$result$subpop_1$mat_unadj_samplevar, test$result$subpop_2$mat_unadj_samplevar)
      mat_unadj_samplevar <- rbind(mat_unadj_samplevar, tmp_unadj_samplevar)
    }
  }
  
  # delete the first line (initializer)
  mat_ltmle <- mat_ltmle[-1, ]
  mat_unadj <- mat_unadj[-1, ]
  mat_ltmle_a1 <- mat_ltmle_a1[-1, ]
  mat_unadj_a1 <- mat_unadj_a1[-1, ]
  mat_ltmle_a0 <- mat_ltmle_a0[-1, ]
  mat_unadj_a0 <- mat_unadj_a0[-1, ]
  mat_ltmle_ICvar <- mat_ltmle_ICvar[-1, ]
  mat_unadj_samplevar <- mat_unadj_samplevar[-1, ]
  
  # swap columns
  if (wait_for_pipeline) {
    curr_id <- c(as.vector(outer(c("Tc_", "tTc_"), 1:K, FUN = paste0)),
                 as.vector(outer(c("T1_", "tT1_"), 1:K, FUN = paste0)),
                 as.vector(outer(c("T2_", "tT2_"), 1:K, FUN = paste0)))
    swap_id <- as.vector(outer(c("Tc_", "tTc_", "T1_", "tT1_", "T2_", "tT2_"), 1:K, FUN = paste0))
  } else {
    curr_id <- as.vector(t(outer(c("Tc_", "T1_", "T2_"), 1:K, FUN = paste0)))
    swap_id <- as.vector(outer(c("Tc_", "T1_", "T2_"), 1:K, FUN = paste0))
  }
  
  colnames(mat_ltmle) <- colnames(mat_unadj) <-
    colnames(mat_ltmle_a1) <- colnames(mat_unadj_a1) <-
    colnames(mat_ltmle_a0) <- colnames(mat_unadj_a0) <-
    colnames(mat_ltmle_ICvar) <- colnames(mat_unadj_samplevar) <- curr_id
  
  mat_ltmle <- mat_ltmle[, swap_id]
  mat_unadj <- mat_unadj[, swap_id]
  mat_ltmle_a1 <- mat_ltmle_a1[, swap_id]
  mat_unadj_a1 <- mat_unadj_a1[, swap_id]
  mat_ltmle_a0 <- mat_ltmle_a0[, swap_id]
  mat_unadj_a0 <- mat_unadj_a0[, swap_id]
  mat_ltmle_ICvar <- mat_ltmle_ICvar[, swap_id]
  mat_unadj_samplevar <- mat_unadj_samplevar[, swap_id]
  
  # Compute mean and covariance of ltmle and unadj
  mean_ltmle <- apply(mat_ltmle, 2, mean)
  cov_ltmle <- cov(mat_ltmle)
  mean_unadj <- apply(mat_unadj, 2, mean)
  cov_unadj <- cov(mat_unadj)
  
  mean_ltmle_a1 <- apply(mat_ltmle_a1, 2, mean)
  cov_ltmle_a1 <- cov(mat_ltmle_a1)
  mean_unadj_a1 <- apply(mat_unadj_a1, 2, mean)
  cov_unadj_a1 <- cov(mat_unadj_a1)
  
  mean_ltmle_a0 <- apply(mat_ltmle_a0, 2, mean)
  cov_ltmle_a0 <- cov(mat_ltmle_a0)
  mean_unadj_a0 <- apply(mat_unadj_a0, 2, mean)
  cov_unadj_a0 <- cov(mat_unadj_a0)
  
  # Standardize ltmle and unadj (divide by standard error, make it Z-statisitc)
  mat_ltmle_std <- mat_ltmle %*% diag(diag(cov_ltmle) ^ (-1/2))
  mat_unadj_std <- mat_unadj %*% diag(diag(cov_unadj) ^ (-1/2))
  
  mat_ltmle_a1_std <- mat_ltmle_a1 %*% diag(diag(cov_ltmle_a1) ^ (-1/2))
  mat_unadj_a1_std <- mat_unadj_a1 %*% diag(diag(cov_unadj_a1) ^ (-1/2))
  
  mat_ltmle_a0_std <- mat_ltmle_a0 %*% diag(diag(cov_ltmle_a0) ^ (-1/2))
  mat_unadj_a0_std <- mat_unadj_a0 %*% diag(diag(cov_unadj_a0) ^ (-1/2))
  
  # Compute mean and covariance of the standardized ones
  mean_ltmle_std <- apply(mat_ltmle_std, 2, mean)
  cov_ltmle_std <- cov(mat_ltmle_std)
  mean_unadj_std <- apply(mat_unadj_std, 2, mean)
  cov_unadj_std <- cov(mat_unadj_std)
  
  mean_ltmle_a1_std <- apply(mat_ltmle_a1_std, 2, mean)
  cov_ltmle_a1_std <- cov(mat_ltmle_a1_std)
  mean_unadj_a1_std <- apply(mat_unadj_a1_std, 2, mean)
  cov_unadj_a1_std <- cov(mat_unadj_a1_std)
  
  mean_ltmle_a0_std <- apply(mat_ltmle_a0_std, 2, mean)
  cov_ltmle_a0_std <- cov(mat_ltmle_a0_std)
  mean_unadj_a0_std <- apply(mat_unadj_a0_std, 2, mean)
  cov_unadj_a0_std <- cov(mat_unadj_a0_std)
  
  # gather the cov/corr results
  result <- list(mat_ltmle = mat_ltmle, mat_unadj = mat_unadj,
                 mean_ltmle = mean_ltmle, mean_unadj = mean_unadj,
                 cov_ltmle = cov_ltmle, cov_unadj = cov_unadj,
                 mat_ltmle_std = mat_ltmle_std, mat_unadj_std = mat_unadj_std,
                 mean_ltmle_std = mean_ltmle_std, mean_unadj_std = mean_unadj_std,
                 cov_ltmle_std = cov_ltmle_std, cov_unadj_std = cov_unadj_std,
                 mat_ltmle_ICvar = mat_ltmle_ICvar,
                 mat_unadj_samplevar = mat_unadj_samplevar,
                 
                 mat_ltmle_a1 = mat_ltmle_a1, mat_unadj_a1 = mat_unadj_a1,
                 mean_ltmle_a1 = mean_ltmle_a1, mean_unadj_a1 = mean_unadj_a1,
                 cov_ltmle_a1 = cov_ltmle_a1, cov_unadj_a1 = cov_unadj_a1,
                 mat_ltmle_a1_std = mat_ltmle_a1_std, mat_unadj_a1_std = mat_unadj_a1_std,
                 mean_ltmle_a1_std = mean_ltmle_a1_std, mean_unadj_a1_std = mean_unadj_a1_std,
                 cov_ltmle_a1_std = cov_ltmle_a1_std, cov_unadj_a1_std = cov_unadj_a1_std,
                 
                 mat_ltmle_a0 = mat_ltmle_a0, mat_unadj_a0 = mat_unadj_a0,
                 mean_ltmle_a0 = mean_ltmle_a0, mean_unadj_a0 = mean_unadj_a0,
                 cov_ltmle_a0 = cov_ltmle_a0, cov_unadj_a0 = cov_unadj_a0,
                 mat_ltmle_a0_std = mat_ltmle_a0_std, mat_unadj_a0_std = mat_unadj_a0_std,
                 mean_ltmle_a0_std = mean_ltmle_a0_std, mean_unadj_a0_std = mean_unadj_a0_std,
                 cov_ltmle_a0_std = cov_ltmle_a0_std, cov_unadj_a0_std = cov_unadj_a0_std)
  return(result)
}

# need to add trials_H01H12
wrapup_gather_results <- function(generic_filename, idx, wait_for_pipeline = FALSE){
  # generic_filename doesn't include H01H02, etc.
  # this returns a list of: trials_H01H02, trials_H11H12, trials_H11H02
  
  names00 <- gen_names(paste0(generic_filename, "H01H02"), idx)
  trials_H01H02 <- gather_results_old(names00, wait_for_pipeline)
  names11 <- gen_names(paste0(generic_filename, "H11H12"), idx)
  trials_H11H12 <- gather_results_old(names11, wait_for_pipeline)
  names10 <- gen_names(paste0(generic_filename, "H11H02"), idx)
  trials_H11H02 <- gather_results_old(names10, wait_for_pipeline)
  
  return(list(trials_H01H02 = trials_H01H02,
              trials_H11H12 = trials_H11H12,
              trials_H11H02 = trials_H11H02))
}

# without recording estimator for mean of each arm (used in theoretical paper)
gather_results_old <- function(filenames, wait_for_pipeline = FALSE){
  
  n <- length(filenames)
  K <- 5
  num_subpop <- 2
  
  ncol <- ifelse(wait_for_pipeline, 2*K*(num_subpop+1), K*(num_subpop+1))
  
  # initialize overall matrix to collect all parallel simulated esitmators
  mat_ltmle <- matrix(NA, nrow = 1, ncol = ncol)
  mat_unadj <- matrix(NA, nrow = 1, ncol = ncol)
  
  # rbind each resulting matrix to the overall matrix  
  for (i in 1:n){
    load(filenames[i])
    tmp_ltmle <- cbind(test$subpop_c$mat_ltmle, test$subpop_1$mat_ltmle, test$subpop_2$mat_ltmle)
    mat_ltmle <- rbind(mat_ltmle, tmp_ltmle)
    tmp_unadj <- cbind(test$subpop_c$mat_unadj, test$subpop_1$mat_unadj, test$subpop_2$mat_unadj)
    mat_unadj <- rbind(mat_unadj, tmp_unadj)
  }
  
  # delete the first line (initializer)
  mat_ltmle <- mat_ltmle[-1, ]
  mat_unadj <- mat_unadj[-1, ]
  
  # swap columns
  if (wait_for_pipeline) {
    curr_id <- as.vector(t(outer(c("Tc_", "tTc_", "T1_", "tT1_", "T2_", "tT2_"), 1:K, FUN = paste0)))
    swap_id <- as.vector(outer(c("Tc_", "tTc_", "T1_", "tT1_", "T2_", "tT2_"), 1:K, FUN = paste0))
  } else {
    curr_id <- as.vector(t(outer(c("Tc_", "T1_", "T2_"), 1:K, FUN = paste0)))
    swap_id <- as.vector(outer(c("Tc_", "T1_", "T2_"), 1:K, FUN = paste0))
  }
  
  colnames(mat_ltmle) <- colnames(mat_unadj) <- curr_id
  
  mat_ltmle <- mat_ltmle[, swap_id]
  mat_unadj <- mat_unadj[, swap_id]
  
  # Compute mean and covariance of ltmle and unadj
  mean_ltmle <- apply(mat_ltmle, 2, mean)
  cov_ltmle <- cov(mat_ltmle)
  mean_unadj <- apply(mat_unadj, 2, mean)
  cov_unadj <- cov(mat_unadj)
  
  # Standardize ltmle and unadj (divide by standard error, make it Z-statisitc)
  mat_ltmle_std <- mat_ltmle %*% diag(diag(cov_ltmle) ^ (-1/2))
  mat_unadj_std <- mat_unadj %*% diag(diag(cov_unadj) ^ (-1/2))
  
  # Compute mean and covariance of the standardized ones
  mean_ltmle_std <- apply(mat_ltmle_std, 2, mean)
  cov_ltmle_std <- cov(mat_ltmle_std)
  mean_unadj_std <- apply(mat_unadj_std, 2, mean)
  cov_unadj_std <- cov(mat_unadj_std)
  
  # gather the cov/corr results
  result <- list(mat_ltmle = mat_ltmle, mat_unadj = mat_unadj,
                 mean_ltmle = mean_ltmle, mean_unadj = mean_unadj,
                 cov_ltmle = cov_ltmle, cov_unadj = cov_unadj,
                 mat_ltmle_std = mat_ltmle_std, mat_unadj_std = mat_unadj_std,
                 mean_ltmle_std = mean_ltmle_std, mean_unadj_std = mean_unadj_std,
                 cov_ltmle_std = cov_ltmle_std, cov_unadj_std = cov_unadj_std)
  return(result)
}
