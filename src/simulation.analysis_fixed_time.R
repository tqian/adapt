### Simulate trials and compute test statistics ###
# Analysis at fixed sample size

## Functions: ##

# simu_trial_fix_ss: simulate nsim trials with analysis at fixed calendar time (specified by timer)
#                    It computes test statistics for H00, H01, H02
# test_stat: compute test statistics (ltmle, unadj) at given time point
# weighted_sum_lists: to compute test statistics for H00, based on H01 and H02 (weighted sum)



library(ltmle)
library(boot)

source("attach_time_to_dataset.R")
source("test_statistic_and_variance.R")

simu_trial_fix_time <- function(nsim, W_to_use, L_to_use = NULL, # L_to_use = NULL means no short-term outcome
                                timer, # time to conduct interim analylsis
                                last_timer_stop_enrollment = FALSE,
                                nmax1, nmax2, # nmax for subpop 1 and 2
                                prevalence_sub1 = 0.469, # prevalence of subpopulation 1 (should be equal to nmax1 / (nmax1 + nmax2))
                                erate_sub1, # erate sub2 will be proportional according to nmax1 and nmax2
                                delay_WL, delay_LY,
                                scenario = c("H11H12", "H01H02", "H11H02", "H01H12"), # if left unspecified, won't modify A
                                dgm, # function of data generating mechanism
                                print_sim_progress = TRUE, # print progress every 10 simulations
                                inter_AW_ltmle = FALSE,
                                enroll_method = c("continuous", "Poisson-homo"), # method to generate enrollment time
                                estimator = c("both", "ltmle", "unadj"),
                                wait_for_pipeline = FALSE,
                                ... # other arguments to be passed into dgm
){
  
  estimator <- match.arg(estimator)
  scenario <- match.arg(scenario)
  
  timer <- as.numeric(timer)
  
  if (last_timer_stop_enrollment) {
    t_stop_enroll <- timer[length(timer)]
    timer[length(timer)] <- timer[length(timer)] + delay_WL + delay_LY + .Machine$double.eps^0.25
  }

  K <- length(timer)
  
  if (wait_for_pipeline) {
    tmp_mat <- matrix(as.numeric(NA), nrow = nsim, ncol = 2 * K)
  } else {
    tmp_mat <- matrix(as.numeric(NA), nrow = nsim, ncol = K)
  }
  
  # matrix to record test statistics
  mat_ltmle_c <- mat_unadj_c <- mat_ltmle_1 <- mat_unadj_1 <- mat_ltmle_2 <- mat_unadj_2 <- tmp_mat  # difference between two arms
  mat_ltmle_c_a1 <- mat_unadj_c_a1 <- mat_ltmle_1_a1 <- mat_unadj_1_a1 <- mat_ltmle_2_a1 <- mat_unadj_2_a1 <- tmp_mat  # arm A = 1
  mat_ltmle_c_a0 <- mat_unadj_c_a0 <- mat_ltmle_1_a0 <- mat_unadj_1_a0 <- mat_ltmle_2_a0 <- mat_unadj_2_a0 <- tmp_mat  # arm A = 0
  mat_ltmle_c_ICvar <- mat_ltmle_1_ICvar <- mat_ltmle_2_ICvar <- tmp_mat # influence curve based variance for treatment effect estimator
  mat_unadj_c_samplevar <- mat_unadj_1_samplevar <- mat_unadj_2_samplevar <- tmp_mat
  
  # matrix to record number enrolled, short-term observed and final observed
  ss_enrolled_c <- ss_L_observed_c <- ss_Y_observed_c <- matrix(as.numeric(NA), nrow = nsim, ncol = K)
  ss_enrolled_1 <- ss_L_observed_1 <- ss_Y_observed_1 <- matrix(as.numeric(NA), nrow = nsim, ncol = K)
  ss_enrolled_2 <- ss_L_observed_2 <- ss_Y_observed_2 <- matrix(as.numeric(NA), nrow = nsim, ncol = K)
  
  # simulate trials to get test statistics
  for (isim in 1:nsim){
    
    if (print_sim_progress == TRUE){
      if( isim %% 10 == 0){
        cat(paste0("\nTrial number ", isim))
      }
    }
    
    dt <- dgm(cumss_sub1 = nmax1, cumss_sub2 = nmax2, ...)
    
    if (scenario == "H01H02") {
      dt$A <- rbinom(nrow(dt), 1, 0.5)
    } else if (scenario == "H11H02") {
      dt$A[dt$subpop == 2] <- rbinom(sum(dt$subpop == 2), 1, 0.5)
    } else if (scenario == "H01H12") {
      dt$A[dt$subpop == 1] <- rbinom(sum(dt$subpop == 1), 1, 0.5)
    }
    
    tmp <- attach_time(dt, nmax1, nmax2, erate_sub1, delay_WL, delay_LY, method = enroll_method)
    dt <- tmp$dt

    if (last_timer_stop_enrollment) {
      dt <- subset(dt, T_enroll <= t_stop_enroll) # remove those who has enrollment time too late
    }

    interim_time <- timer
    
    dt_combined <- dt
    dt_subpop1 <- subset(dt, subpop == 1)
    dt_subpop2 <- subset(dt, subpop == 2)
    
    for (k in 1:K) {
      time <- as.numeric(interim_time[k])
      ss_enrolled_c[isim, k] <- sum(dt_combined$T_enroll <= time)
      ss_L_observed_c[isim, k] <- sum(dt_combined$T_shortterm <= time)
      ss_Y_observed_c[isim, k] <- sum(dt_combined$T_final <= time)
      ss_enrolled_1[isim, k] <- sum(dt_subpop1$T_enroll <= time)
      ss_L_observed_1[isim, k] <- sum(dt_subpop1$T_shortterm <= time)
      ss_Y_observed_1[isim, k] <- sum(dt_subpop1$T_final <= time)
      ss_enrolled_2[isim, k] <- sum(dt_subpop2$T_enroll <= time)
      ss_L_observed_2[isim, k] <- sum(dt_subpop2$T_shortterm <= time)
      ss_Y_observed_2[isim, k] <- sum(dt_subpop2$T_final <= time)
    }
    
    
    ests_1 <- test_stat(dt_subpop1, interim_time, W_to_use, L_to_use,
                        inter_AW_ltmle = inter_AW_ltmle, estimator = estimator, wait_for_pipeline = wait_for_pipeline)
    ests_2 <- test_stat(dt_subpop2, interim_time, W_to_use, L_to_use,
                        inter_AW_ltmle = inter_AW_ltmle, estimator = estimator, wait_for_pipeline = wait_for_pipeline)
    
    weight_sub1 <- prevalence_sub1
    ests_c <- weighted_sum_lists(l1 = ests_1, l2 = ests_2, w1 = weight_sub1, w2 = 1- weight_sub1)
    
    mat_ltmle_c[isim, ] <- ests_c$ests["ltmle", ]
    mat_unadj_c[isim, ] <- ests_c$ests["unadj", ]
    mat_ltmle_1[isim, ] <- ests_1$ests["ltmle", ]
    mat_unadj_1[isim, ] <- ests_1$ests["unadj", ]
    mat_ltmle_2[isim, ] <- ests_2$ests["ltmle", ]
    mat_unadj_2[isim, ] <- ests_2$ests["unadj", ]
    
    mat_ltmle_c_ICvar[isim, ] <- ests_c$ests["ltmle_ICvar", ]
    mat_ltmle_1_ICvar[isim, ] <- ests_1$ests["ltmle_ICvar", ]
    mat_ltmle_2_ICvar[isim, ] <- ests_2$ests["ltmle_ICvar", ]
    
    mat_unadj_c_samplevar[isim, ] <- ests_c$ests["unadj_samplevar", ]
    mat_unadj_1_samplevar[isim, ] <- ests_1$ests["unadj_samplevar", ]
    mat_unadj_2_samplevar[isim, ] <- ests_2$ests["unadj_samplevar", ]
    
    mat_ltmle_c_a1[isim, ] <- ests_c$ests_a1["ltmle", ]
    mat_unadj_c_a1[isim, ] <- ests_c$ests_a1["unadj", ]
    mat_ltmle_1_a1[isim, ] <- ests_1$ests_a1["ltmle", ]
    mat_unadj_1_a1[isim, ] <- ests_1$ests_a1["unadj", ]
    mat_ltmle_2_a1[isim, ] <- ests_2$ests_a1["ltmle", ]
    mat_unadj_2_a1[isim, ] <- ests_2$ests_a1["unadj", ]
    
    mat_ltmle_c_a0[isim, ] <- ests_c$ests_a0["ltmle", ]
    mat_unadj_c_a0[isim, ] <- ests_c$ests_a0["unadj", ]
    mat_ltmle_1_a0[isim, ] <- ests_1$ests_a0["ltmle", ]
    mat_unadj_1_a0[isim, ] <- ests_1$ests_a0["unadj", ]
    mat_ltmle_2_a0[isim, ] <- ests_2$ests_a0["ltmle", ]
    mat_unadj_2_a0[isim, ] <- ests_2$ests_a0["unadj", ]
  }
  
  result <- list(subpop_c = list(), subpop_1 = list(), subpop_2 = list())  
  for (i in 1:3){
    if (i == 1) {
      mat_ltmle <- mat_ltmle_c
      mat_unadj <- mat_unadj_c
      mat_ltmle_ICvar <- mat_ltmle_c_ICvar
      mat_unadj_samplevar <- mat_unadj_c_samplevar
    } else if (i == 2) {
      mat_ltmle <- mat_ltmle_1
      mat_unadj <- mat_unadj_1
      mat_ltmle_ICvar <- mat_ltmle_1_ICvar
      mat_unadj_samplevar <- mat_unadj_1_samplevar
    } else if (i == 3) {
      mat_ltmle <- mat_ltmle_2
      mat_unadj <- mat_unadj_2
      mat_ltmle_ICvar <- mat_ltmle_2_ICvar
      mat_unadj_samplevar <- mat_unadj_2_samplevar
    }    
    result[[i]] <- list(mat_ltmle = mat_ltmle, mat_unadj = mat_unadj,
                        mat_ltmle_ICvar = mat_ltmle_ICvar, mat_unadj_samplevar = mat_unadj_samplevar)
  }
  
  result_a1 <- list(subpop_c = list(), subpop_1 = list(), subpop_2 = list())  
  for (i in 1:3){
    if (i == 1) {
      mat_ltmle <- mat_ltmle_c_a1
      mat_unadj <- mat_unadj_c_a1
    } else if (i == 2) {
      mat_ltmle <- mat_ltmle_1_a1
      mat_unadj <- mat_unadj_1_a1
    } else if (i == 3) {
      mat_ltmle <- mat_ltmle_2_a1
      mat_unadj <- mat_unadj_2_a1
    }    
    result_a1[[i]] <- list(mat_ltmle = mat_ltmle, mat_unadj = mat_unadj)
  }
  
  result_a0 <- list(subpop_c = list(), subpop_1 = list(), subpop_2 = list())  
  for (i in 1:3){
    if (i == 1) {
      mat_ltmle <- mat_ltmle_c_a0
      mat_unadj <- mat_unadj_c_a0
    } else if (i == 2) {
      mat_ltmle <- mat_ltmle_1_a0
      mat_unadj <- mat_unadj_1_a0
    } else if (i == 3) {
      mat_ltmle <- mat_ltmle_2_a0
      mat_unadj <- mat_unadj_2_a0
    }    
    result_a0[[i]] <- list(mat_ltmle = mat_ltmle, mat_unadj = mat_unadj)
  }
  
  return(list(result = result, result_a1 = result_a1, result_a0 = result_a0,
              ss_enrolled_c = ss_enrolled_c, ss_L_observed_c = ss_L_observed_c, ss_Y_observed_c = ss_Y_observed_c,
              ss_enrolled_1 = ss_enrolled_1, ss_L_observed_1 = ss_L_observed_1, ss_Y_observed_1 = ss_Y_observed_1,
              ss_enrolled_2 = ss_enrolled_2, ss_L_observed_2 = ss_L_observed_2, ss_Y_observed_2 = ss_Y_observed_2))
}


weighted_sum_lists <- function(l1, l2, w1 = 0.5, w2 = 0.5, var_term = c("ltmle_ICvar", "unadj_samplevar")){
  l = list()
  for(i in 1:length(l1)){
    l[[i]] <- w1 * l1[[i]] + w2 * l2[[i]]
    for (variable in var_term) {
      l[[i]][variable, ] <- w1^2 * l1[[i]][variable, ] + w2^2 * l2[[i]][variable, ]
    }
  }
  names(l) <- names(l1)
  return(l)
}