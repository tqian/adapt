### Evaluate simulated trials when testing for H00, H01 and H02, using error spending approach ###
# Analysis at fixed time

## Functions: ##

# eval_trial_fix_time_H00H01H02: takes in test statistics, error to spend at each stage;
#                       computes testing boundaries, and outputs the trial result (which hypo rejected, and when rejected)

# stopping_detail_H00H01H02: decides which and when each hypo is rejected, and when each enrollment is stopped


library(mvtnorm)

##### change this line if different directory structure #####
source("boundary_error_spending.R")


eval_trial_fix_time_H00H01H02 <- function(trials_H01H02, trials_H11H12, trials_H11H02, trials_H01H12,
                                          true_trt_eff_H1 = NULL, # for calculating bias - currently not being used
                                          nmax1, nmax2, timer, 
                                          erate_sub1 = 70, delay_WL = 30/365, delay_LY = 150/365,
                                          alphas_user = NULL, betas_user = NULL, # errors at each stage
                                          alpha_c = NULL, alpha_1 = NULL, alpha_2,
                                          f_err = NULL, # error spending functions
                                          binding_fut = FALSE,
                                          modify_l = "specify-vector",
                                          l_1k = NULL, l_2k = NULL, # two vectors for l_{1,k}, l_{2,k}
                                          estimator = c("both", "ltmle", "unadj"),
                                          print = c("basic", "more", "no")){
  
  
  estimator <- match.arg(estimator)
  
  ### decide print level ###
  
  print_level <- match.arg(print)
  
  ### some constants ###
  nsim <- nrow(trials_H01H02$mat_ltmle)
  K <- length(timer)
  stg_always_stop_2 <- K
  
  ### For calculating ESS and duration ###
  analysis_time <- as.numeric(timer)
  erate_sub2 <- erate_sub1 * nmax2 / nmax1
  
  ### compute relative efficiency ###
  
  if (estimator == "both") {
    
    ## RE for difference in mean
    RE_H01H02 <- diag(trials_H01H02$cov_unadj) / diag(trials_H01H02$cov_ltmle)
    RE_table_H01H02 <- matrix(RE_H01H02, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_table_H01H02) <- c("combined", "subpop1", "subpop2")
    colnames(RE_table_H01H02) <- paste("stg", 1:K)
    
    RE_H11H12 <- diag(trials_H11H12$cov_unadj) / diag(trials_H11H12$cov_ltmle)
    RE_table_H11H12 <- matrix(RE_H11H12, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_table_H11H12) <- c("combined", "subpop1", "subpop2")
    colnames(RE_table_H11H12) <- paste("stg", 1:K)
    
    RE_H11H02 <- diag(trials_H11H02$cov_unadj) / diag(trials_H11H02$cov_ltmle)
    RE_table_H11H02 <- matrix(RE_H11H02, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_table_H11H02) <- c("combined", "subpop1", "subpop2")
    colnames(RE_table_H11H02) <- paste("stg", 1:K)
    
    RE_H01H12 <- diag(trials_H01H12$cov_unadj) / diag(trials_H01H12$cov_ltmle)
    RE_table_H01H12 <- matrix(RE_H01H12, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_table_H01H12) <- c("combined", "subpop1", "subpop2")
    colnames(RE_table_H01H12) <- paste("stg", 1:K)
    
    ## RE for mean in A = 1
    RE_a1_H01H02 <- diag(trials_H01H02$cov_unadj_a1) / diag(trials_H01H02$cov_ltmle_a1)
    RE_a1_table_H01H02 <- matrix(RE_a1_H01H02, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a1_table_H01H02) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a1_table_H01H02) <- paste("stg", 1:K)
    
    RE_a1_H11H12 <- diag(trials_H11H12$cov_unadj_a1) / diag(trials_H11H12$cov_ltmle_a1)
    RE_a1_table_H11H12 <- matrix(RE_a1_H11H12, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a1_table_H11H12) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a1_table_H11H12) <- paste("stg", 1:K)
    
    RE_a1_H11H02 <- diag(trials_H11H02$cov_unadj_a1) / diag(trials_H11H02$cov_ltmle_a1)
    RE_a1_table_H11H02 <- matrix(RE_a1_H11H02, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a1_table_H11H02) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a1_table_H11H02) <- paste("stg", 1:K)
    
    RE_a1_H01H12 <- diag(trials_H01H12$cov_unadj_a1) / diag(trials_H01H12$cov_ltmle_a1)
    RE_a1_table_H01H12 <- matrix(RE_a1_H01H12, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a1_table_H01H12) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a1_table_H01H12) <- paste("stg", 1:K)
    
    ## RE for mean in A = 0
    RE_a0_H01H02 <- diag(trials_H01H02$cov_unadj_a0) / diag(trials_H01H02$cov_ltmle_a0)
    RE_a0_table_H01H02 <- matrix(RE_a0_H01H02, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a0_table_H01H02) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a0_table_H01H02) <- paste("stg", 1:K)
    
    RE_a0_H11H12 <- diag(trials_H11H12$cov_unadj_a0) / diag(trials_H11H12$cov_ltmle_a0)
    RE_a0_table_H11H12 <- matrix(RE_a0_H11H12, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a0_table_H11H12) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a0_table_H11H12) <- paste("stg", 1:K)
    
    RE_a0_H11H02 <- diag(trials_H11H02$cov_unadj_a0) / diag(trials_H11H02$cov_ltmle_a0)
    RE_a0_table_H11H02 <- matrix(RE_a0_H11H02, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a0_table_H11H02) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a0_table_H11H02) <- paste("stg", 1:K)
    
    RE_a0_H01H12 <- diag(trials_H01H12$cov_unadj_a0) / diag(trials_H01H12$cov_ltmle_a0)
    RE_a0_table_H01H12 <- matrix(RE_a0_H01H12, nrow = 3, ncol = K, byrow = FALSE)
    rownames(RE_a0_table_H01H12) <- c("combined", "subpop1", "subpop2")
    colnames(RE_a0_table_H01H12) <- paste("stg", 1:K)    
    
  } # end if (estimator == "both")
  
  if (estimator %in% c("both", "ltmle")) {
    
    ###############
    #### ltmle ####
    ###############
    
    idx_c <- seq(from = 1, by = 3, length = K)
    idx_1 <- seq(from = 2, by = 3, length = K)
    idx_2 <- seq(from = 3, by = 3, length = K)
    
    ### get error spending boundary ###
    alphas_c_ltmle <- errs_alpha_only(var_vector = diag(trials_H01H02$cov_ltmle[idx_c, idx_c]),
                                      alpha = alpha_c, alphas = alphas_user[1:K],
                                      f_err = f_err, 
                                      always_stop_at = stg_always_stop_2)
    alphas_1_ltmle <- errs_alpha_only(var_vector = diag(trials_H01H02$cov_ltmle[idx_1, idx_1]),
                                      alpha = alpha_1, alphas = alphas_user[(K+1):(2*K)],
                                      f_err = f_err)
    alphas_2_ltmle <- errs_alpha_only(var_vector = diag(trials_H01H02$cov_ltmle[idx_2, idx_2]),
                                      alpha = alpha_2, alphas = alphas_user[(2*K+1):(3*K)],
                                      f_err = f_err, 
                                      always_stop_at = stg_always_stop_2)
    
    # 4 sets of boundaries, only difference is the covariance matrix (should be very similar)
    bdry_ltmle_H01H02 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_ltmle_std,
                                               trials_H01H02$cov_ltmle_std,
                                               alphas_c_ltmle, alphas_1_ltmle, alphas_2_ltmle,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)    
    bdry_ltmle_H11H02 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_ltmle_std,
                                               trials_H11H02$cov_ltmle_std,
                                               alphas_c_ltmle, alphas_1_ltmle, alphas_2_ltmle,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)    
    bdry_ltmle_H01H12 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_ltmle_std,
                                               trials_H01H12$cov_ltmle_std,
                                               alphas_c_ltmle, alphas_1_ltmle, alphas_2_ltmle,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)
    bdry_ltmle_H11H12 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_ltmle_std,
                                               trials_H11H12$cov_ltmle_std,
                                               alphas_c_ltmle, alphas_1_ltmle, alphas_2_ltmle,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)
    

    
    if (print_level == "more"){
      cat("ltmle, alphas spent at each stage:\n")
      print(round(alphas_c_ltmle,5))
      print(round(alphas_1_ltmle,5))
      print(round(alphas_2_ltmle,5))
      
      cat("ltmle, error spending boundaries (H01H02):\n")
      print(lapply(bdry_ltmle_H01H02, round, digits = 2))
      cat("ltmle, error spending boundaries (H11H12):\n")
      print(lapply(bdry_ltmle_H11H12, round, digits = 2))
      cat("ltmle, error spending boundaries (H11H02):\n")
      print(lapply(bdry_ltmle_H11H02, round, digits = 2))
      cat("ltmle, error spending boundaries (H01H12):\n")
      print(lapply(bdry_ltmle_H01H12, round, digits = 2))
    }
    
    ### get trial results ###
    
    results_ltmle_H01H02_b <- stopping_detail_H00H01H02(trials_H01H02$mat_ltmle_std, bdry_ltmle_H01H02, binding_fut = T)
    results_ltmle_H11H12_b <- stopping_detail_H00H01H02(trials_H11H12$mat_ltmle_std, bdry_ltmle_H11H12, binding_fut = T)
    results_ltmle_H11H02_b <- stopping_detail_H00H01H02(trials_H11H02$mat_ltmle_std, bdry_ltmle_H11H02, binding_fut = T)
    results_ltmle_H01H12_b <- stopping_detail_H00H01H02(trials_H01H12$mat_ltmle_std, bdry_ltmle_H01H12, binding_fut = T)
    
    results_ltmle_H01H02_nb <- stopping_detail_H00H01H02(trials_H01H02$mat_ltmle_std, bdry_ltmle_H01H02, binding_fut = F)
    results_ltmle_H11H12_nb <- stopping_detail_H00H01H02(trials_H11H12$mat_ltmle_std, bdry_ltmle_H11H12, binding_fut = F)
    results_ltmle_H11H02_nb <- stopping_detail_H00H01H02(trials_H11H02$mat_ltmle_std, bdry_ltmle_H11H02, binding_fut = F)
    results_ltmle_H01H12_nb <- stopping_detail_H00H01H02(trials_H01H12$mat_ltmle_std, bdry_ltmle_H01H12, binding_fut = F)
    
    ### Type I error and power ###
    
    # For type I error, use nonbinding results (conservative)
    # For power, use binding results (conservative)
    typeIerror_ltmle_H01H02 <- sum(results_ltmle_H01H02_nb$reject_H01 + results_ltmle_H01H02_nb$reject_H02 + results_ltmle_H01H02_nb$reject_H0C > 0) / nsim
    typeIerror_ltmle_H11H02 <- sum(results_ltmle_H11H02_nb$reject_H02 > 0) / nsim
    typeIerror_ltmle_H01H12 <- sum(results_ltmle_H01H12_nb$reject_H01 > 0) / nsim
    
    # power of rejecting H0C when H11H12
    power_ltmle_c <- sum(results_ltmle_H11H12_b$reject_H0C) / nsim
    # power of rejecting H01 when H11H02
    power_ltmle_1 <- sum(results_ltmle_H11H02_b$reject_H01) / nsim
    # power of rejecting H02 when H01H12
    power_ltmle_2 <- sum(results_ltmle_H01H12_b$reject_H02) / nsim
    
    scns <- c("H01H02", "H01H12", "H11H02", "H11H12")
    
    for (bindingtype in c('b', 'nb')) {
      tmp_reject_prob <- matrix(nrow = 4, ncol = 7)
      rownames(tmp_reject_prob) <- scns
      colnames(tmp_reject_prob) <- c("Pow_H0C", "Pow_H01", "Pow_H02", "Pow_H01_and_H0C", "Pow_H02_and_H0C", "Pow_all", "Pow_any")
      
      for (scn in scns) {
        eval(parse(text = paste0("tmp_results_ltmle <- results_ltmle_", scn, "_", bindingtype)))
        tmp_reject_prob[scn, "Pow_H0C"] <- sum(tmp_results_ltmle$reject_H0C)
        tmp_reject_prob[scn, "Pow_H01"] <- sum(tmp_results_ltmle$reject_H01)
        tmp_reject_prob[scn, "Pow_H02"] <- sum(tmp_results_ltmle$reject_H02)
        tmp_reject_prob[scn, "Pow_H01_and_H0C"] <- sum(tmp_results_ltmle$reject_H01 + tmp_results_ltmle$reject_H0C == 2)
        tmp_reject_prob[scn, "Pow_H02_and_H0C"] <- sum(tmp_results_ltmle$reject_H02 + tmp_results_ltmle$reject_H0C == 2)
        tmp_reject_prob[scn, "Pow_all"] <- sum(tmp_results_ltmle$reject_H0C + tmp_results_ltmle$reject_H01 + tmp_results_ltmle$reject_H02 == 3)
        tmp_reject_prob[scn, "Pow_any"] <- sum(tmp_results_ltmle$reject_H0C + tmp_results_ltmle$reject_H01 + tmp_results_ltmle$reject_H02 > 0)
      }
      tmp_reject_prob <- tmp_reject_prob / nsim
      eval(parse(text = paste0("reject_prob_ltmle_", bindingtype, " <- tmp_reject_prob")))
    }

    
    
    ### expected sample size (ESS), average duration, and percentage stopping at each stage ###  
    
    computeESS_timer <- function(simulated_trials, result_list){
      
      rowid <- 1:nsim
      
      colid <- result_list$subpop_1_stopped_when
      ind <- (colid - 1)* nsim + rowid
      ESS_1 <- mean(simulated_trials$ss_enrolled_1[ind])
      
      colid <- result_list$subpop_2_stopped_when
      ind <- (colid - 1)* nsim + rowid
      ESS_2 <- mean(simulated_trials$ss_enrolled_2[ind])
      
      ESS_c <- ESS_1 + ESS_2
      
      ESS <- matrix(c(ESS_c, ESS_1, ESS_2), nrow = 3)
      return(ESS)
    }
    
    computeDuration <- function(result_list){
      duration_1 <- mean(analysis_time[result_list$subpop_1_stopped_when])
      duration_2 <- mean(analysis_time[result_list$subpop_2_stopped_when])
      duration_all <- max(duration_1, duration_2)
      
      duration <- matrix(c(duration_all, duration_1, duration_2), nrow = 3)
      return(duration)
    }
    
    stop_table <- function(result_list){
      tbl <- matrix(-999, ncol = K + 1, nrow = 6)
      colnames(tbl) <- c(paste("stg", 1:K), "total")
      rownames(tbl) <- c("reject H0C", "accept H0C", "reject H01", "accept H01", "reject H02", "accept H02")
      
      for (i in 1:K) {
        tbl["reject H0C", i] <- mean(min(result_list$subpop_1_stopped_when, result_list$subpop_2_stopped_when) == i & result_list$reject_H0C == 1)
        tbl["accept H0C", i] <- mean(min(result_list$subpop_1_stopped_when, result_list$subpop_2_stopped_when) == i & result_list$reject_H0C == 0)
        tbl["reject H01", i] <- mean(result_list$subpop_1_stopped_when == i & result_list$reject_H01 == 1)
        tbl["accept H01", i] <- mean(result_list$subpop_1_stopped_when == i & result_list$reject_H01 == 0)
        tbl["reject H02", i] <- mean(result_list$subpop_2_stopped_when == i & result_list$reject_H02 == 1)
        tbl["accept H02", i] <- mean(result_list$subpop_2_stopped_when == i & result_list$reject_H02 == 0)
      }
      for (j in 1:nrow(tbl)) tbl[j, K+1] <- sum(tbl[j, 1:K])
      return(tbl)
    }
    
    results_all_ltmle <- list(results_ltmle_H01H02_b, results_ltmle_H11H12_b, results_ltmle_H11H02_b, results_ltmle_H01H12_b,
                              results_ltmle_H01H02_nb, results_ltmle_H11H12_nb, results_ltmle_H11H02_nb, results_ltmle_H01H12_nb)
    
    ESS_ltmle <- matrix(as.numeric(NA), nrow = 3, ncol = 8)
    rownames(ESS_ltmle) <- c("ESS_c", "ESS_1", "ESS_2")
    colnames(ESS_ltmle) <- c("H01H02_b", "H11H12_b", "H11H02_b", "H01H12_b",
                             "H01H02_nb", "H11H12_nb", "H11H02_nb",  "H01H12_nb")
    for (iscn in 1:8){
      if (iscn %in% c(1,5)) simulated_trials <- trials_H01H02
      if (iscn %in% c(2,6)) simulated_trials <- trials_H11H12
      if (iscn %in% c(3,7)) simulated_trials <- trials_H11H02
      if (iscn %in% c(4,8)) simulated_trials <- trials_H01H12
      ESS_ltmle[, iscn] <- computeESS_timer(simulated_trials, results_all_ltmle[[iscn]])
    }
    
    duration_ltmle <- sapply(results_all_ltmle, computeDuration)
    rownames(duration_ltmle) <- c("duration_all", "duration_1", "duration_2")
    colnames(duration_ltmle) <- c("H01H02_b", "H11H12_b", "H11H02_b", "H01H12_b",
                                  "H01H02_nb", "H11H12_nb", "H11H02_nb",  "H01H12_nb")
    
    stop_tables_ltmle <- lapply(results_all_ltmle, stop_table)
    names(stop_tables_ltmle) <- c("H01H02_b", "H11H12_b", "H11H02_b", "H01H12_b",
                                  "H01H02_nb", "H11H12_nb", "H11H02_nb",  "H01H12_nb")
    
  } # end if (estimator %in% c("both", "ltmle"))
  
  
  if (estimator %in% c("both", "unadj")) {
    
    ###############
    #### unadj ####
    ###############
    
    idx_c <- seq(from = 1, by = 3, length = K)
    idx_1 <- seq(from = 2, by = 3, length = K)
    idx_2 <- seq(from = 3, by = 3, length = K)
    
    ### get error spending boundary ###
    alphas_c_unadj <- errs_alpha_only(var_vector = diag(trials_H01H02$cov_unadj[idx_c, idx_c]),
                                      alpha = alpha_c, alphas = alphas_user[1:K],
                                      f_err = f_err, 
                                      always_stop_at = stg_always_stop_2)
    alphas_1_unadj <- errs_alpha_only(var_vector = diag(trials_H01H02$cov_unadj[idx_1, idx_1]),
                                      alpha = alpha_1, alphas = alphas_user[(K+1):(2*K)],
                                      f_err = f_err)
    alphas_2_unadj <- errs_alpha_only(var_vector = diag(trials_H01H02$cov_unadj[idx_2, idx_2]),
                                      alpha = alpha_2, alphas = alphas_user[(2*K+1):(3*K)],
                                      f_err = f_err, 
                                      always_stop_at = stg_always_stop_2)
    
    # 4 sets of boundaries, only difference is the covariance matrix (should be very similar)
    bdry_unadj_H01H02 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_unadj_std,
                                               trials_H01H02$cov_unadj_std,
                                               alphas_c_unadj, alphas_1_unadj, alphas_2_unadj,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)    
    bdry_unadj_H11H02 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_unadj_std,
                                               trials_H11H02$cov_unadj_std,
                                               alphas_c_unadj, alphas_1_unadj, alphas_2_unadj,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)    
    bdry_unadj_H01H12 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_unadj_std,
                                               trials_H01H12$cov_unadj_std,
                                               alphas_c_unadj, alphas_1_unadj, alphas_2_unadj,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)
    bdry_unadj_H11H12 <- err_sp_bdry_H00H01H02(trials_H01H02$mean_unadj_std,
                                               trials_H11H12$cov_unadj_std,
                                               alphas_c_unadj, alphas_1_unadj, alphas_2_unadj,
                                               modify_l = modify_l,
                                               l_1k = l_1k, l_2k = l_2k)
    
    
    
    if (print_level == "more"){
      cat("unadj, alphas spent at each stage:\n")
      print(round(alphas_c_unadj,5))
      print(round(alphas_1_unadj,5))
      print(round(alphas_2_unadj,5))
      
      cat("unadj, error spending boundaries (H01H02):\n")
      print(lapply(bdry_unadj_H01H02, round, digits = 2))
      cat("unadj, error spending boundaries (H11H12):\n")
      print(lapply(bdry_unadj_H11H12, round, digits = 2))
      cat("unadj, error spending boundaries (H11H02):\n")
      print(lapply(bdry_unadj_H11H02, round, digits = 2))
      cat("unadj, error spending boundaries (H01H12):\n")
      print(lapply(bdry_unadj_H01H12, round, digits = 2))
    }
    
    ### get trial results ###
    
    results_unadj_H01H02_b <- stopping_detail_H00H01H02(trials_H01H02$mat_unadj_std, bdry_unadj_H01H02, binding_fut = T)
    results_unadj_H11H12_b <- stopping_detail_H00H01H02(trials_H11H12$mat_unadj_std, bdry_unadj_H11H12, binding_fut = T)
    results_unadj_H11H02_b <- stopping_detail_H00H01H02(trials_H11H02$mat_unadj_std, bdry_unadj_H11H02, binding_fut = T)
    results_unadj_H01H12_b <- stopping_detail_H00H01H02(trials_H01H12$mat_unadj_std, bdry_unadj_H01H12, binding_fut = T)
    
    results_unadj_H01H02_nb <- stopping_detail_H00H01H02(trials_H01H02$mat_unadj_std, bdry_unadj_H01H02, binding_fut = F)
    results_unadj_H11H12_nb <- stopping_detail_H00H01H02(trials_H11H12$mat_unadj_std, bdry_unadj_H11H12, binding_fut = F)
    results_unadj_H11H02_nb <- stopping_detail_H00H01H02(trials_H11H02$mat_unadj_std, bdry_unadj_H11H02, binding_fut = F)
    results_unadj_H01H12_nb <- stopping_detail_H00H01H02(trials_H01H12$mat_unadj_std, bdry_unadj_H01H12, binding_fut = F)
    
    ### Type I error and power ###
    
    # For type I error, use nonbinding results (conservative)
    # For power, use binding results (conservative)
    typeIerror_unadj_H01H02 <- sum(results_unadj_H01H02_nb$reject_H01 + results_unadj_H01H02_nb$reject_H02 + results_unadj_H01H02_nb$reject_H0C > 0) / nsim
    typeIerror_unadj_H11H02 <- sum(results_unadj_H11H02_nb$reject_H02 > 0) / nsim
    typeIerror_unadj_H01H12 <- sum(results_unadj_H01H12_nb$reject_H01 > 0) / nsim
    
    # power of rejecting H0C when H11H12
    power_unadj_c <- sum(results_unadj_H11H12_b$reject_H0C) / nsim
    # power of rejecting H01 when H11H02
    power_unadj_1 <- sum(results_unadj_H11H02_b$reject_H01) / nsim
    # power of rejecting H02 when H01H12
    power_unadj_2 <- sum(results_unadj_H01H12_b$reject_H02) / nsim
    
    scns <- c("H01H02", "H01H12", "H11H02", "H11H12")
    
    for (bindingtype in c('b', 'nb')) {
      tmp_reject_prob <- matrix(nrow = 4, ncol = 7)
      rownames(tmp_reject_prob) <- scns
      colnames(tmp_reject_prob) <- c("Pow_H0C", "Pow_H01", "Pow_H02", "Pow_H01_and_H0C", "Pow_H02_and_H0C", "Pow_all", "Pow_any")
      
      for (scn in scns) {
        eval(parse(text = paste0("tmp_results_unadj <- results_unadj_", scn, "_", bindingtype)))
        tmp_reject_prob[scn, "Pow_H0C"] <- sum(tmp_results_unadj$reject_H0C)
        tmp_reject_prob[scn, "Pow_H01"] <- sum(tmp_results_unadj$reject_H01)
        tmp_reject_prob[scn, "Pow_H02"] <- sum(tmp_results_unadj$reject_H02)
        tmp_reject_prob[scn, "Pow_H01_and_H0C"] <- sum(tmp_results_unadj$reject_H01 + tmp_results_unadj$reject_H0C == 2)
        tmp_reject_prob[scn, "Pow_H02_and_H0C"] <- sum(tmp_results_unadj$reject_H02 + tmp_results_unadj$reject_H0C == 2)
        tmp_reject_prob[scn, "Pow_all"] <- sum(tmp_results_unadj$reject_H0C + tmp_results_unadj$reject_H01 + tmp_results_unadj$reject_H02 == 3)
        tmp_reject_prob[scn, "Pow_any"] <- sum(tmp_results_unadj$reject_H0C + tmp_results_unadj$reject_H01 + tmp_results_unadj$reject_H02 > 0)
      }
      tmp_reject_prob <- tmp_reject_prob / nsim
      eval(parse(text = paste0("reject_prob_unadj_", bindingtype, " <- tmp_reject_prob")))
    }    
    
    ### expected sample size (ESS), average duration, and percentage stopping at each stage ###  
    
    computeESS_timer <- function(simulated_trials, result_list){
      
      rowid <- 1:nsim
      
      colid <- result_list$subpop_1_stopped_when
      ind <- (colid - 1)* nsim + rowid
      ESS_1 <- mean(simulated_trials$ss_enrolled_1[ind])
      
      colid <- result_list$subpop_2_stopped_when
      ind <- (colid - 1)* nsim + rowid
      ESS_2 <- mean(simulated_trials$ss_enrolled_2[ind])
      
      ESS_c <- ESS_1 + ESS_2
      
      ESS <- matrix(c(ESS_c, ESS_1, ESS_2), nrow = 3)
      return(ESS)
    }
    
    computeDuration <- function(result_list){
      duration_1 <- mean(analysis_time[result_list$subpop_1_stopped_when])
      duration_2 <- mean(analysis_time[result_list$subpop_2_stopped_when])
      duration_all <- max(duration_1, duration_2)
      
      duration <- matrix(c(duration_all, duration_1, duration_2), nrow = 3)
      return(duration)
    }
    
    stop_table <- function(result_list){
      tbl <- matrix(-999, ncol = K + 1, nrow = 6)
      colnames(tbl) <- c(paste("stg", 1:K), "total")
      rownames(tbl) <- c("reject H0C", "accept H0C", "reject H01", "accept H01", "reject H02", "accept H02")
      
      for (i in 1:K) {
        tbl["reject H0C", i] <- mean(result_list$subpop_c_rejected_when == i & result_list$reject_H0C == 1)
        tbl["accept H0C", i] <- 0 # doesn't mean anything
        tbl["reject H01", i] <- mean(result_list$subpop_1_stopped_when == i & result_list$reject_H01 == 1)
        tbl["accept H01", i] <- mean(result_list$subpop_1_stopped_when == i & result_list$reject_H01 == 0)
        tbl["reject H02", i] <- mean(result_list$subpop_2_stopped_when == i & result_list$reject_H02 == 1)
        tbl["accept H02", i] <- mean(result_list$subpop_2_stopped_when == i & result_list$reject_H02 == 0)
      }
      for (j in 1:nrow(tbl)) tbl[j, K+1] <- sum(tbl[j, 1:K])
      return(tbl)
    }
    
    results_all_unadj <- list(results_unadj_H01H02_b, results_unadj_H11H12_b, results_unadj_H11H02_b, results_unadj_H01H12_b,
                              results_unadj_H01H02_nb, results_unadj_H11H12_nb, results_unadj_H11H02_nb, results_unadj_H01H12_nb)
    
    ESS_unadj <- matrix(as.numeric(NA), nrow = 3, ncol = 8)
    rownames(ESS_unadj) <- c("ESS_c", "ESS_1", "ESS_2")
    colnames(ESS_unadj) <- c("H01H02_b", "H11H12_b", "H11H02_b", "H01H12_b",
                             "H01H02_nb", "H11H12_nb", "H11H02_nb",  "H01H12_nb")
    for (iscn in 1:8){
      if (iscn %in% c(1,5)) simulated_trials <- trials_H01H02
      if (iscn %in% c(2,6)) simulated_trials <- trials_H11H12
      if (iscn %in% c(3,7)) simulated_trials <- trials_H11H02
      if (iscn %in% c(4,8)) simulated_trials <- trials_H01H12
      ESS_unadj[, iscn] <- computeESS_timer(simulated_trials, results_all_unadj[[iscn]])
    }
    
    duration_unadj <- sapply(results_all_unadj, computeDuration)
    rownames(duration_unadj) <- c("duration_all", "duration_1", "duration_2")
    colnames(duration_unadj) <- c("H01H02_b", "H11H12_b", "H11H02_b", "H01H12_b",
                                  "H01H02_nb", "H11H12_nb", "H11H02_nb",  "H01H12_nb")
    
    stop_tables_unadj <- lapply(results_all_unadj, stop_table)
    names(stop_tables_unadj) <- c("H01H02_b", "H11H12_b", "H11H02_b", "H01H12_b",
                                  "H01H02_nb", "H11H12_nb", "H11H02_nb",  "H01H12_nb")
    
  } # end if (estimator %in% c("both", "unadj"))

  
  
  ### Output ###
  if (estimator == "both"){
    
    output <- list(estimator = estimator,
                   nmax1 = nmax1, nmax2 = nmax2, timer = timer,
                   nsim = nsim,
                   typeIerror = matrix(c(typeIerror_ltmle_H01H02, typeIerror_ltmle_H11H02, typeIerror_ltmle_H01H12,
                                         typeIerror_unadj_H01H02, typeIerror_unadj_H11H02, typeIerror_unadj_H01H12), nrow = 2, ncol = 3, byrow = TRUE),
                   power = matrix(c(power_ltmle_c, power_ltmle_1, power_ltmle_2,
                             power_unadj_c, power_unadj_1, power_unadj_2), nrow = 2, ncol = 3, byrow = TRUE),
                   reject_prob_ltmle_b = reject_prob_ltmle_b, reject_prob_ltmle_nb = reject_prob_ltmle_nb,
                   reject_prob_unadj_b = reject_prob_unadj_b, reject_prob_unadj_nb = reject_prob_unadj_nb,
                   RE = list(RE_table = list(RE_table_H01H02 = RE_table_H01H02,
                                             RE_table_H11H12 = RE_table_H11H12,
                                             RE_table_H11H02 = RE_table_H11H02,
                                             RE_table_H01H12 = RE_table_H01H12),
                             RE_a1_table = list(RE_a1_table_H01H02 = RE_a1_table_H01H02,
                                                RE_a1_table_H11H12 = RE_a1_table_H11H12,
                                                RE_a1_table_H11H02 = RE_a1_table_H11H02,
                                                RE_a1_table_H01H12 = RE_a1_table_H01H12),
                             RE_a0_table = list(RE_a0_table_H01H02 = RE_a0_table_H01H02,
                                                RE_a0_table_H11H12 = RE_a0_table_H11H12,
                                                RE_a0_table_H11H02 = RE_a0_table_H11H02,
                                                RE_a0_table_H01H12 = RE_a0_table_H01H12)),                             
                   expected_sample_size = list(ltmle = ESS_ltmle, unadj = ESS_unadj),
                   average_duration = list(ltmle = duration_ltmle, unadj = duration_unadj),
                   stop_tables_ltmle = stop_tables_ltmle,
                   stop_tables_unadj = stop_tables_unadj,
                   alphas_c_ltmle = alphas_c_ltmle,
                   alphas_1_ltmle = alphas_1_ltmle,
                   alphas_c_unadj = alphas_c_unadj,
                   alphas_1_unadj = alphas_1_unadj,
                   bdry_ltmle = list(H01H02 = bdry_ltmle_H01H02, H11H12 = bdry_ltmle_H11H12, H11H02 = bdry_ltmle_H11H02, H01H12 = bdry_ltmle_H01H12),
                   bdry_unadj = list(H01H02 = bdry_unadj_H01H02, H11H12 = bdry_unadj_H11H12, H11H02 = bdry_unadj_H11H02, H01H12 = bdry_unadj_H01H12),
                   mean_ltmle = list(H01H02 = trials_H01H02$mean_ltmle,
                                     H01H12 = trials_H01H12$mean_ltmle,
                                     H11H02 = trials_H11H02$mean_ltmle,
                                     H11H12 = trials_H11H12$mean_ltmle),
                   mean_unadj = list(H01H02 = trials_H01H02$mean_unadj,
                                     H01H12 = trials_H01H12$mean_unadj,
                                     H11H02 = trials_H11H02$mean_unadj,
                                     H11H12 = trials_H11H12$mean_unadj),
                   mean_ltmle_std = list(H01H02 = trials_H01H02$mean_ltmle_std,
                                         H01H12 = trials_H01H12$mean_ltmle_std,
                                         H11H02 = trials_H11H02$mean_ltmle_std,
                                         H11H12 = trials_H11H12$mean_ltmle_std),
                   mean_unadj_std = list(H01H02 = trials_H01H02$mean_unadj_std,
                                         H01H12 = trials_H01H12$mean_unadj_std,
                                         H11H02 = trials_H11H02$mean_unadj_std,
                                         H11H12 = trials_H11H12$mean_unadj_std),
                   cov_ltmle = list(H01H02 = trials_H01H02$cov_ltmle,
                                    H01H12 = trials_H01H12$cov_ltmle,
                                    H11H02 = trials_H11H02$cov_ltmle,
                                    H11H12 = trials_H11H12$cov_ltmle),
                   cov_unadj = list(H01H02 = trials_H01H02$cov_unadj,
                                    H01H12 = trials_H01H12$cov_unadj,
                                    H11H02 = trials_H11H02$cov_unadj,
                                    H11H12 = trials_H11H12$cov_unadj),
                   cov_ltmle_std = list(H01H02 = trials_H01H02$cov_ltmle_std,
                                        H01H12 = trials_H01H12$cov_ltmle_std,
                                        H11H02 = trials_H11H02$cov_ltmle_std,
                                        H11H12 = trials_H11H12$cov_ltmle_std),
                   cov_unadj_std = list(H01H02 = trials_H01H02$cov_unadj_std,
                                        H01H12 = trials_H01H12$cov_unadj_std,
                                        H11H02 = trials_H11H02$cov_unadj_std,
                                        H11H12 = trials_H11H12$cov_unadj_std))
    colnames(output$typeIerror) <- c("H01H02", "H11H02", "H01H12")
    rownames(output$typeIerror) <- c("ltmle", "unadj")
    colnames(output$power) <- c("H0C", "H01", "H02")
    rownames(output$power) <- c("ltmle", "unadj")
    
  } else if (estimator == "unadj") {
    
    output <- list(estimator = estimator,
                   nmax1 = nmax1, nmax2 = nmax2, timer = timer,
                   nsim = nsim,
                   typeIerror = c(typeIerror_unadj_H01H02, typeIerror_unadj_H11H02, typeIerror_unadj_H01H12),
                   power = c(power_unadj_c, power_unadj_1, power_unadj_2),
                   reject_prob_unadj_b = reject_prob_unadj_b, reject_prob_unadj_nb = reject_prob_unadj_nb,
                   expected_sample_size = list(unadj = ESS_unadj),
                   average_duration = list(unadj = duration_unadj),
                   stop_tables_unadj = stop_tables_unadj,
                   alphas_c_unadj = alphas_c_unadj,
                   alphas_1_unadj = alphas_1_unadj,
                   bdry_unadj = list(H01H02 = bdry_unadj_H01H02, H11H12 = bdry_unadj_H11H12, H11H02 = bdry_unadj_H11H02, H01H12 = bdry_unadj_H01H12),
                   mean_unadj = list(H01H02 = trials_H01H02$mean_unadj,
                                     H01H12 = trials_H01H12$mean_unadj,
                                     H11H02 = trials_H11H02$mean_unadj,
                                     H11H12 = trials_H11H12$mean_unadj),
                   mean_unadj_std = list(H01H02 = trials_H01H02$mean_unadj_std,
                                         H01H12 = trials_H01H12$mean_unadj_std,
                                         H11H02 = trials_H11H02$mean_unadj_std,
                                         H11H12 = trials_H11H12$mean_unadj_std),
                   cov_unadj = list(H01H02 = trials_H01H02$cov_unadj,
                                    H01H12 = trials_H01H12$cov_unadj,
                                    H11H02 = trials_H11H02$cov_unadj,
                                    H11H12 = trials_H11H12$cov_unadj),
                   cov_unadj_std = list(H01H02 = trials_H01H02$cov_unadj_std,
                                        H01H12 = trials_H01H12$cov_unadj_std,
                                        H11H02 = trials_H11H02$cov_unadj_std,
                                        H11H12 = trials_H11H12$cov_unadj_std))
    names(output$typeIerror) <- c("H01H02", "H11H02", "H01H12")
    names(output$power) <- c("H0C", "H01", "H02")
  } else if (estimator == "ltmle") {
    
    output <- list(estimator = estimator,
                   nmax1 = nmax1, nmax2 = nmax2, timer = timer,
                   nsim = nsim,
                   typeIerror = c(typeIerror_ltmle_H01H02, typeIerror_ltmle_H11H02, typeIerror_ltmle_H01H12),
                   power = c(power_ltmle_c, power_ltmle_1, power_ltmle_2),
                   reject_prob_ltmle_b = reject_prob_ltmle_b, reject_prob_ltmle_nb = reject_prob_ltmle_nb,
                   expected_sample_size = list(ltmle = ESS_ltmle),
                   average_duration = list(ltmle = duration_ltmle),
                   stop_tables_ltmle = stop_tables_ltmle,
                   alphas_c_ltmle = alphas_c_ltmle,
                   alphas_1_ltmle = alphas_1_ltmle,
                   bdry_ltmle = list(H01H02 = bdry_ltmle_H01H02, H11H12 = bdry_ltmle_H11H12, H11H02 = bdry_ltmle_H11H02, H01H12 = bdry_ltmle_H01H12),
                   mean_ltmle = list(H01H02 = trials_H01H02$mean_ltmle,
                                     H01H12 = trials_H01H12$mean_ltmle,
                                     H11H02 = trials_H11H02$mean_ltmle,
                                     H11H12 = trials_H11H12$mean_ltmle),
                   mean_ltmle_std = list(H01H02 = trials_H01H02$mean_ltmle_std,
                                         H01H12 = trials_H01H12$mean_ltmle_std,
                                         H11H02 = trials_H11H02$mean_ltmle_std,
                                         H11H12 = trials_H11H12$mean_ltmle_std),
                   cov_ltmle = list(H01H02 = trials_H01H02$cov_ltmle,
                                    H01H12 = trials_H01H12$cov_ltmle,
                                    H11H02 = trials_H11H02$cov_ltmle,
                                    H11H12 = trials_H11H12$cov_ltmle),
                   cov_ltmle_std = list(H01H02 = trials_H01H02$cov_ltmle_std,
                                        H01H12 = trials_H01H12$cov_ltmle_std,
                                        H11H02 = trials_H11H02$cov_ltmle_std,
                                        H11H12 = trials_H11H12$cov_ltmle_std))
    names(output$typeIerror) <- c("H01H02", "H11H02", "H01H12")
    names(output$power) <- c("H0C", "H01", "H02")
  }
  
  return(output)
}



stopping_detail_H00H01H02 <- function(mat_zstats, bdry, stg_always_stop_2 = NULL, binding_fut = TRUE){
  
  K <- ncol(mat_zstats) / 3
  
  idx_c <- seq(from = 1, by = 3, length = K)
  idx_1 <- seq(from = 2, by = 3, length = K)
  idx_2 <- seq(from = 3, by = 3, length = K)
  
  Z_c <- mat_zstats[, idx_c]
  Z_1 <- mat_zstats[, idx_1]
  Z_2 <- mat_zstats[, idx_2]
  
  u_c <- bdry$u_c
  u_1 <- bdry$u_1
  u_2 <- bdry$u_2
  
  l_1 <- bdry$l_1
  l_2 <- bdry$l_2
  
  reject_H0C <- reject_H01 <- reject_H02 <- rep(0, nrow(mat_zstats))
  subpop_1_stopped <- subpop_2_stopped <- rep(0, nrow(mat_zstats))
  subpop_c_rejected_when <- subpop_1_stopped_when <- subpop_2_stopped_when <- rep(NA, nrow(mat_zstats))
  
  if (is.null(stg_always_stop_2)){
    stg_always_stop_2 <- K
  }
  
  for(stage in 1:K)
  {
    
    # Determine if any new events where a null hypothesis is rejected for efficacy:
    
    reject_H01 <- ifelse((!subpop_1_stopped) & Z_1[, stage] > u_1[stage], 1, reject_H01)
    
    reject_H02 <- ifelse((!subpop_2_stopped) & Z_2[, stage] > u_2[stage], 1 ,reject_H02)
    
    reject_H0C <- ifelse( (reject_H01 & reject_H02) |
                            ((!subpop_1_stopped) & (!subpop_2_stopped) & Z_c[, stage] > u_c[stage]),
                          1, reject_H0C)
    
    if (binding_fut) {
      subpop_1_stopped <- ifelse(stage == K | reject_H01 | (Z_1[, stage] < l_1[stage]), 1, subpop_1_stopped)    
      subpop_2_stopped <- ifelse(stage == stg_always_stop_2 | reject_H02 | (Z_2[, stage] < l_2[stage]), 1, subpop_2_stopped)
    } else {
      subpop_1_stopped <- ifelse(stage == K | reject_H01, 1, subpop_1_stopped)    
      subpop_2_stopped <- ifelse(stage == stg_always_stop_2 | reject_H02, 1, subpop_2_stopped)
    }
    
    subpop_c_rejected_when <- ifelse(is.na(subpop_c_rejected_when) & reject_H0C, stage, subpop_c_rejected_when)
    subpop_1_stopped_when <- ifelse(is.na(subpop_1_stopped_when) & subpop_1_stopped, stage, subpop_1_stopped_when)
    subpop_2_stopped_when <- ifelse(is.na(subpop_2_stopped_when) & subpop_2_stopped, stage, subpop_2_stopped_when)
    
  }
  
  return(list(reject_H01 = reject_H01, reject_H02 = reject_H02, reject_H0C = reject_H0C,
              subpop_c_rejected_when = subpop_c_rejected_when,
              subpop_1_stopped_when = subpop_1_stopped_when, subpop_2_stopped_when = subpop_2_stopped_when))
}
