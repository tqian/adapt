
# last row of beta.W, beta.WL, betaL.W are standard deviations of the error term in linear regression

# Format of the coefficients: (each should have 2 columns, for subpop 1 and subpop 2)
#
# beta.W:  coef for W on Y when there is no L
#          each column has intercept, coef for W
# beta.WL: coef for W,L on Y
#          each column has intercept, coef for W, and coef for L (last row)
# betaL.W: coef for W on Y
#          each column has intercept, coef for W
# betaA:  coef for A on Y
# betaL.A: coef for A on L
#
# sigmaY: standard deviation of error for Y
# sigmaL: standard deviation of error for L




dgm <- function(cumss_sub1, cumss_sub2,
                betas = NULL, # coefficients beta_W and beta_L
                alphas = NULL, # coefficients alpha_W
                betaA = c(0.2353647, 0.1765766), # coefficient of A for Y and L
                alphaA = 0.2121, # obtained in AD_dgm_design_new.R
                sigmaY, sigmaL = NULL, # standard deviation of the error terms in regression model
                trt_eff = 0.4242,
                data, S = "subpop", W = c("W1", "W2", "W3", "W4", "W5"),
                AD_dgm_scn = c("H01H02", "H11H02", "H01H12", "H11H12"),
                exW = FALSE, # set W exogenous after generating Y and L
                exL = FALSE, # set L exogenous after generating Y and L
                p_resetL = 0,
                rand.p = 0.5 # randomization probability (of treatment)
                ) {
  

  if (is.null(trt_eff)) {
    trt_eff <- 0.4242
  }
  
  if (length(p_resetL) == 1) {
    p_resetL <- rep(p_resetL, 2)
  }
  
  meanL <- c(mean(data$L[data[, S] == 1]), mean(data$L[data[, S] == 2]))
  sdL <- c(sd(data$L[data[, S] == 1]), sd(data$L[data[, S] == 2]))
  
  # simulation scenario (whether setting A exogenous)
  scn <- match.arg(AD_dgm_scn)
  
  ### Calculate the maximum sample size
  n.sub <- c(cumss_sub1[length(cumss_sub1)], cumss_sub2[length(cumss_sub2)])
  
  ### Create numeric value for the ID
  data$ID <- seq(1,nrow(data))
  
  ### Generate the data separately for each subgroup
  for(i in c(1,2)) { # i: subpopulation 1, 2
    
    ### Sample W with replacement within subset S
    S.ids <- sample(data$ID[data[, S] == i], size = n.sub[i], replace = T)
    
    ### Generate dataframe for each subset and assign treatment with prob rand.p
    data.new <- data[S.ids,]
    data.new$A <- rbinom(n.sub[i], 1, prob = rand.p)
    
    no_effect <- (i == 1 & scn %in% c("H01H02", "H01H12")) | (i == 2 & scn %in% c("H01H02", "H11H02"))
    
    if (no_effect) {
      m.new <- as.matrix(cbind(rep(1, nrow(data.new)), data.new[, W]))
      data.new$L <- m.new %*% alphas[, i] + rnorm(nrow(data.new), mean = 0, sd = sigmaL[i])
      m.new.L <- as.matrix(cbind(rep(1, nrow(data.new)), data.new[, W], data.new$L))
      data.new$Y <- m.new.L %*% betas[, i] + rnorm(nrow(data.new), mean = 0, sd = sigmaY[i])
    } else { # A has effect
      m.new <- as.matrix(cbind(rep(1,nrow(data.new)), data.new[,W], data.new$A))
      data.new$L <- m.new %*% c(alphas[, i], alphaA) + rnorm(nrow(data.new), mean = 0, sd = sigmaL[i])
      m.new.L <- as.matrix(cbind(rep(1,nrow(data.new)), data.new[,W], data.new$L, data.new$A))
      data.new$Y <- m.new.L %*% c(betas[, i], betaA[i]) + rnorm(nrow(data.new), mean = 0, sd = sigmaY[i])
    }
    
    if (exW) {
      for (iW in 1:length(W)) {
        data.new[, W[iW]] <- rnorm(length(data.new[, W[iW]]), mean = 0, sd = sd(data.new[, W[iW]]))
      }
    }
  
    
    if (exL) {
      data.new[, 'L'] <- rnorm(length(data.new[, 'L']), mean = meanL[i], sd = sdL[i])
    }
    
    if (p_resetL[i] > 0) {
      id_reset <- sample.int(nrow(data.new), round(p_resetL[i] * nrow(data.new)))
      data.new[id_reset, 'L'] <- rnorm(length(data.new[id_reset, 'L']), mean = meanL[i], sd = sdL[i])
    }

    if(i==1) df_1 <- data.new
    if(i==2) df_2 <- data.new
  }
  df <- rbind(df_1,df_2)
  return(df)
}


### a wrapper for simu_eval_trial() ###

AD_dgm_wrapper <- function(setup) {
  dt <- dgm(cumss_sub1 = setup$cumss_sub1,
            cumss_sub2 = setup$cumss_sub2,
            beta.W = setup$beta.W,
            beta.WL = setup$beta.WL,
            betaL.W = setup$betaL.W,
            betaA = setup$betaA,
            betaL.A = setup$betaL.A,
            sigmaY = setup$sigmaY,
            sigmaL = setup$sigmaL,
            trt_eff = setup$.trt_eff,
            data = setup$data,
            exW = setup$exW,
            exL = setup$exL,
            AD_dgm_scn = setup$effset,
            includeL = setup$includeL)
  
#   if (setup$effset %in% c("H01H02", "H11H02")) {
#     dt$A[dt$subpop == 2] <- rbinom(sum(dt$subpop == 2), 1, 0.5)
#   }
#   if (setup$effset %in% c("H01H02", "H01H12")) {
#     dt$A[dt$subpop == 1] <- rbinom(sum(dt$subpop == 1), 1, 0.5)
#   }

  return(dt)
}
