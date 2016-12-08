### Compute test statistic (both, ltmle, unadj) and the variance at given time point(s) ###

## Functions: ##

# test_stat: computes test statistics (treatment effect, mean of each arm) and variance (ICvar or bootstrapped for ltmle) at interim_time
# weighted_sum_lists: to compute test statistics for H00, based on H01 and H02 (weighted sum)


# To work on:
# 1) test_stat_pkgltmle
#   1.1) wait_for_pipeline & inter_AW_ltmle case
#   1.2) bootstrap to estimate variance of test statistics

# 2) test_stat_ownltmle
#   2.1) inter_AW_ltmle case
#   2.2) wait_for_pipeline case
#   2.3) bootstrap to estimate variance of test statistics

test_stat <- function(dt, interim_time, W_to_use, L_to_use, wait_for_pipeline = FALSE,
                               inter_AW_ltmle = FALSE, estimator = c("both", "ltmle", "unadj"),
                               ltmle_method = c("package", "own_code")
                               #                       ,nboot = nbootstrap_var
){
  
  estimator <- match.arg(estimator)
  ltmle_method <- match.arg(ltmle_method)
  K <- length(interim_time)
  
  if (wait_for_pipeline) {
    tmp <- matrix(NA, nrow = 4, ncol = 2*K)
    colnames(tmp) <- as.vector(outer(c("T_", "tT_"), 1:K, paste0))
    rownames(tmp) <- c("ltmle", "unadj", "ltmle_ICvar", "unadj_samplevar") # influence curve based variance estimate for ltmle treatment effect
    ests <- ests_a1 <- ests_a0 <- tmp
  } else {
    tmp <- matrix(NA, nrow = 4, ncol = K)
    colnames(tmp) <- paste0("T_", 1:K)
    rownames(tmp) <- c("ltmle", "unadj", "ltmle_ICvar", "unadj_samplevar")
    ests <- ests_a1 <- ests_a0 <- tmp
  }
  
  for (istage in 1:K){
    ## Interim analysis
    
    idx <- ifelse(wait_for_pipeline, 2*istage - 1, istage)
    
    time <- as.numeric(interim_time[istage])
    C0 <- as.numeric(dt$T_enroll <= time)
    C1 <- as.numeric(dt$T_shortterm <= time)
    C2 <- as.numeric(dt$T_final <= time)
    
    # construct the data frame to be passed into ltmle()
    dt_ltmle <- data.frame(C0, dt[, W_to_use], dt$A,
                           C1, dt[, L_to_use],
                           C2, dt$Y)
    colnames(dt_ltmle) <- c("C0", W_to_use, "A",
                            "C1", L_to_use,
                            "C2", "Y")
    
    if (inter_AW_ltmle) {
      # Currently this only works for:
      # interaction terms are W3*A and W4*A, W_to_use = c("W3", "W4"), L_to_use = NULL
      gform <- c("C0 ~ 1", "A ~ W3 + W4", "C1 ~ W3*A + W4*A", "C2 ~ W3*A + W4*A")
      Qform <- c(W3="Q.kplus1 ~ 1", W4="Q.kplus1 ~ W3", 
                 #                  L="Q.kplus1 ~ W3*A + W4*A", 
                 Y="Q.kplus1 ~ W3*A + W4*A")
    } else {
      gform <- NULL
      Qform <- NULL
    }
    
    
    if (estimator %in% c("both", "ltmle")) {      
      # here W_to_use are to the left of Anode, so they will be treated as baseline variables
      
      if (ltmle_method == "package") {
        suppressMessages(ltmle.fit <- ltmle(dt_ltmle, 
                                            Anodes = "A",
                                            Cnodes = c("C0","C1","C2"),
                                            Lnodes = c(W_to_use, L_to_use),
                                            Ynodes = "Y",
                                            abar = list(1, 0),
                                            estimate.time = FALSE,
                                            gform = gform, Qform = Qform,
                                            IC.variance.only = TRUE))
        Yrange <- attr(summary(ltmle.fit)$transformOutcome, "Yrange")
        Yleft <- Yrange[1]
        Yrange <- Yrange[2] - Yrange[1]
        
        ests["ltmle", idx] <- as.numeric(summary(ltmle.fit)$effect.measures$ATE$estimate) * Yrange
        ests_a1["ltmle", idx] <- as.numeric(summary(ltmle.fit)$effect.measures$treatment$estimate) * Yrange + Yleft
        ests_a0["ltmle", idx] <- as.numeric(summary(ltmle.fit)$effect.measures$control$estimate) * Yrange + Yleft
        ests["ltmle_ICvar", idx] <- as.numeric(summary(ltmle.fit)$effect.measures$ATE$std.dev)^2 * Yrange^2
        ests_a1["ltmle_ICvar", idx] <- as.numeric(summary(ltmle.fit)$effect.measures$treatment$std.dev)^2 * Yrange^2
        ests_a0["ltmle_ICvar", idx] <- as.numeric(summary(ltmle.fit)$effect.measures$control$std.dev)^2 * Yrange^2
      } else {
        ltmle.fit <- ltmle_ATE(dt_ltmle, Anodes = "A",  Cnodes = c("C0","C1","C2"), Wnodes = W_to_use, Lnodes = L_to_use, Ynodes = "Y")
        ests["ltmle", idx] <- ltmle.fit["ATE"]
        ests_a1["ltmle", idx] <- ltmle.fit["treatment"]
        ests_a0["ltmle", idx] <- ltmle.fit["control"]
      }
      
    }
    
    if (estimator %in% c("both", "unadj")) {
      ests["unadj", idx] <- mean(subset(dt_ltmle, (C2==1) & (A==1))$Y) - mean(subset(dt_ltmle, (C2==1) & (A==0))$Y)
      ests_a1["unadj", idx] <- mean(subset(dt_ltmle, (C2==1) & (A==1))$Y)
      ests_a0["unadj", idx] <- mean(subset(dt_ltmle, (C2==1) & (A==0))$Y)
      n1 <- length(subset(dt_ltmle, (C2==1) & (A==1))$Y)
      var_Y1 <- var(subset(dt_ltmle, (C2==1) & (A==1))$Y)
      n0 <- length(subset(dt_ltmle, (C2==1) & (A==0))$Y)
      var_Y0 <- var(subset(dt_ltmle, (C2==1) & (A==0))$Y)
      ests["unadj_samplevar", idx] <- var_Y1 / n1 + var_Y0 / n0
      ests_a1["unadj_samplevar", idx] <- var_Y1 / n1
      ests_a0["unadj_samplevar", idx] <- var_Y0 / n0
    }
    
    if (wait_for_pipeline){
      # estimator when current pipelines are finished (assuming not enrolling more patients)
      
      # decision analysis should be based on all the already-enrolled people
      dt_tmp <- dt[which(C0 == 1), ]
      
      # construct the data frame to be passed into ltmle()
      # Here don't need short-term outcome
      dt_ltmle <- data.frame(dt_tmp[, W_to_use], dt_tmp$A, dt_tmp$Y)
      colnames(dt_ltmle) <- c(W_to_use, "A", "Y")
      
      # if (inter_AW_ltmle)
      # need to work on this
      
      if (estimator %in% c("both", "ltmle")) {
        suppressMessages(ltmle.fit <- ltmle(dt_ltmle, 
                                            Anodes = "A",
                                            Ynodes = "Y",
                                            abar = list(1, 0),
                                            estimate.time = FALSE,
                                            IC.variance.only = TRUE))
        Yrange <- attr(summary(ltmle.fit)$transformOutcome, "Yrange")
        Yleft <- Yrange[1]
        Yrange <- Yrange[2] - Yrange[1]
        
        ests["ltmle", idx + 1] <- as.numeric(summary(ltmle.fit)$effect.measures$ATE$estimate) * Yrange
        ests["ltmle_ICvar", idx + 1] <- as.numeric(summary(ltmle.fit)$effect.measures$ATE$std.dev)^2 * Yrange^2
        ests_a1["ltmle", idx + 1] <- as.numeric(summary(ltmle.fit)$effect.measures$treatment$estimate) * Yrange + Yleft
        ests_a0["ltmle", idx + 1] <- as.numeric(summary(ltmle.fit)$effect.measures$control$estimate) * Yrange + Yleft
      }
      
      if (estimator %in% c("both", "unadj")) {
        ests["unadj", idx + 1] <- mean(subset(dt_ltmle, A==1)$Y) - mean(subset(dt_ltmle, A==0)$Y)
        ests_a1["unadj", idx + 1] <- mean(subset(dt_ltmle, A==1)$Y)
        ests_a0["unadj", idx + 1] <- mean(subset(dt_ltmle, A==0)$Y)
        n1 <- length(subset(dt_ltmle, A==1)$Y)
        var_Y1 <- var(subset(dt_ltmle, A==1)$Y)
        n0 <- length(subset(dt_ltmle, A==0)$Y)
        var_Y0 <- var(subset(dt_ltmle, A==0)$Y)
        ests["unadj_samplevar", idx + 1] <- var_Y1 / n1 + var_Y0 / n0
        ests_a1["unadj_samplevar", idx + 1] <- var_Y1 / n1
        ests_a0["unadj_samplevar", idx + 1] <- var_Y0 / n0
      }
    }
  }
  
  return(list(ests = ests, ests_a1 = ests_a1, ests_a0 = ests_a0))
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