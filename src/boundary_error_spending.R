### compute efficacy boundaries based on error spending approach ###

## Functions: ##

# err_sp_bdry_H00H01H02.wrapper: wrapper function to be used in simu_eval_trial
# err_sp_bdry_H00H01H02: compute efficacy boundaries using mean, vcov, and alpha allocation at each stage
# err_sp_bdry_H00H01: compute both efficacy and futility boundaries (using approach in Hampson and Jennison 2013 paper)

# errs_alpha_only: compute alpha  allocation at each stage, using variance vector, total alpha, spending function
# errs: compute alpha and beta allocation at each stage, using variance vector, total alpha, total beta, spending functions

# power_ESF: error spending function
# gen_power_ESF: generate error spending function with different parameters
# .f_err, .g_err: ESF with power = 2

# conservative error spending boundaries based on discussion with Michael on 1/20/2015

library(mvtnorm)


### 2015.07.30 ###
# wrapper function, to be used in simu_eval_trial

err_sp_bdry_H00H01H02.wrapper <- function(setup, ests) {
  K <- setup$K
  
  bdry <- err_sp_bdry_H00H01H02(mvmean0 = rep(0, 3*K),
                                mvcov0 = ests$vcov_std,
                                alphas_c = setup$alphas_c,
                                alphas_1 = setup$alphas_1,
                                alphas_2 = setup$alphas_2)
  return(list(u_c = bdry$u_c, u_1 = bdry$u_1, u_2 = bdry$u_2))
}

### End 2015.07.30 ###

# This is for testing H00, H01 and H02, for the simulation paper
err_sp_bdry_H00H01H02 <- function(mvmean0, mvcov0, # mean vector and covariance matrix of standardized test statistics
                               alphas_c, alphas_1, alphas_2, # type I errors of H00, H01, H02
                               modify_l = "specify-vector",
                               l_1k = NULL, l_2k = NULL # two vectors for l_{1,k}, l_{2,k}
){

  alphas <- as.vector(rbind(alphas_c, alphas_1, alphas_2))
  
  K <- length(alphas_c)
  
  u <- as.numeric(rep(NA, length = length(alphas)))
  
  for(k in 1:length(u)){
    # do binary search to get all the u and l's
    
    if (k == 1){ # first stage for H0C and H01
      u[k] <- qnorm(alphas[k], mvmean0[k], sqrt(mvcov0[k,k]), lower.tail = FALSE)      
    } else { # for H0C, not stage 1
      if (alphas[k] > 0) {
        u[k] <- uniroot_trycatch( # equation (12)
          function(x){
            pmvnorm(lower = c(rep(-Inf, k-1), x), upper = c(u[1:(k-1)], Inf),
                    mean = mvmean0[1:k], sigma = mvcov0[1:k, 1:k],
                    algorithm = GenzBretz(abseps = 0.000001, maxpts = 500000)) - alphas[k] },
          interval = c(1, 4), extendInt = "yes", maxiter = 1000)
      } else {
        u[k] <- Inf
      }
      
    }
    
  }
  
  idx_c <- seq(from = 1, by = 3, length = K)
  idx_1 <- seq(from = 2, by = 3, length = K)
  idx_2 <- seq(from = 3, by = 3, length = K)
  
  u_c <- u[idx_c]
  u_1 <- u[idx_1]
  u_2 <- u[idx_2]
  
  l_1 <- l_1k
  l_2 <- l_2k
    
  bdry <- list(u_c = u_c, u_1 = u_1, u_2 = u_2,
               l_1 = l_1, l_2 = l_2)
  
  return(bdry)
}



# This is for testing H00 and H01, for the theoretical paper
err_sp_bdry_H00H01 <- function(mvmean0, mvcov0, mvmean1, mvcov1, # mean vector and covariance matrix of standardized test statistics
                        errs_c, errs_1, # errs of H0C and H01, returned by errs()
                        modify_l = c("all-zero", "conservative-inf", "l2-neg-inf"), # choose the way to modify l
                        binding_fut = FALSE # in calculating u, whether assume binding futility boundaries or not
){
  
  modify_l <- match.arg(modify_l)
  
  alphas_c <- errs_c[1, ]
  betas_c <- errs_c[2, ]
  alphas_1 <- errs_1[1, ]
  betas_1 <- errs_1[2, ]
  
  alphas <- as.vector(rbind(alphas_c, alphas_1))
  betas <- as.vector(rbind(betas_c, betas_1))
  
  K <- length(alphas_c)
  
  # Testing for H0C will stop at this stage
  stg_always_stop_c <- which(alphas_c == 0)[1] - 1
  if (length(stg_always_stop_c) == 0) {
    stg_always_stop_c <- K
  }
  
  u <- l <- as.numeric(rep(NA, length = 2*K))
  
  for(k in 1:(2*K)){
    # do binary search to get all the u and l's
    
    if (k <= 2*1){ # first stage for H0C and H01
      u[k] <- qnorm(alphas[k], mvmean0[k], sqrt(mvcov0[k,k]), lower.tail = FALSE) # equation before (12)
      l[k] <- qnorm(betas[k], mvmean1[k], sqrt(mvcov1[k,k])) # equation before (12)
      
    } else if (k %% 2 == 1){ # for H0C, not stage 1
      
      if (k > 2*stg_always_stop_c) {
        u[k] <- Inf
        l[k] <- -Inf
      } else {
        
        idx <- 1:k # idx for computing u[k] and l[k]        
        if (binding_fut) {
          u[k] <- uniroot_trycatch( # equation (12)
            function(x){
              pmvnorm(lower = c(l[1:(k-1)], x), upper = c(u[1:(k-1)], Inf),
                      mean = mvmean0[idx], sigma = mvcov0[idx, idx],
                      algorithm = GenzBretz(abseps = 0.000001, maxpts = 500000)) - alphas[k] },
            interval = c(-10, 10), extendInt = "yes", maxiter = 200)
        } else {
          u[k] <- uniroot_trycatch( # equation (12)
            function(x){
              pmvnorm(lower = c(rep(-Inf, k-1), x), upper = c(u[1:(k-1)], Inf),
                      mean = mvmean0[idx], sigma = mvcov0[idx, idx],
                      algorithm = GenzBretz(abseps = 0.000001, maxpts = 500000)) - alphas[k] },
            interval = c(-10, 10), extendInt = "yes", maxiter = 200)
        }
        
        if (k %in% c(2*stg_always_stop_c-1, 2*K-1)) {
          l[k] <- Inf
        } else {
          l[k] <- uniroot_trycatch( # equation (13)
            function(x){
              pmvnorm(lower = c(l[1:(k-1)], -Inf), upper = c(u[1:(k-1)], x),
                      mean = mvmean1[idx], sigma = mvcov1[idx, idx],
                      algorithm = GenzBretz(abseps = 0.000001, maxpts = 500000)) - betas[k] },
            interval = c(-10, 10), extendInt = "yes", maxiter = 200)
        }
        
      }
    } else { # for H01, not stage 1
      
      idx <- 1:k
      
      if (binding_fut) {
        u[k] <- uniroot_trycatch( # equation (12)
          function(x){
            pmvnorm(lower = c(l[1:(k-1)], x), upper = c(u[1:(k-1)], Inf),
                    mean = mvmean0[idx], sigma = mvcov0[idx, idx],
                    algorithm = GenzBretz(abseps = 0.000001, maxpts = 500000)) - alphas[k] },
          interval = c(-10, 10), extendInt = "yes", maxiter = 200)
      } else {
        u[k] <- uniroot_trycatch( # equation (12)
          function(x){
            pmvnorm(lower = c(rep(-Inf, k-1), x), upper = c(u[1:(k-1)], Inf),
                    mean = mvmean0[idx], sigma = mvcov0[idx, idx],
                    algorithm = GenzBretz(abseps = 0.000001, maxpts = 500000)) - alphas[k] },
          interval = c(-10, 10), extendInt = "yes", maxiter = 200)
      }
      
      if (k == 2*K) {
        l[k] <- u[k]
      } else {
        l[k] <- uniroot_trycatch( # equation (13)
          function(x){
            pmvnorm(lower = c(l[1:(k-1)], -Inf), upper = c(u[1:(k-1)], x),
                    mean = mvmean1[idx], sigma = mvcov1[idx, idx],
                    algorithm = GenzBretz(abseps = 0.000001, maxpts = 500000)) - betas[k] },
          interval = c(-10, 10), extendInt = "yes", maxiter = 200)
      }
      
    }      
  }
  
  if (modify_l == "conservative-inf") {
    l[is.nan(l)] <- -Inf
  } else if (modify_l == "all-zero") {
    l[-c(5,10)] <- 0
  } else if (modify_l == "l2-neg-inf") {
    l[-c(5,10)] <- 0
    l[c(1,3)] <- -Inf
  }
  
  #   l[is.nan(u)] <- Inf
  #   u[is.nan(u)] <- Inf
  #   l[is.nan(l)] <- -Inf
  #   l[-c(5,10)] <- 0
  #   u[is.infinite(l)] <- Inf # This line might be too conservative
  #                           # This is for the situation where at some stage u is computable but l is not.
  #                           # e.g. nmax = 440
  
  bdry <- list(u = u, l = l)
  
  return(bdry)
}


uniroot_trycatch <- function(f, interval, ...){
  
  solution <- tryCatch(
{
  uniroot(f, interval, ...)$root
},
error = function(cond) {
  message("\nCatched error in uniroot:")
  message(cond)
  message("\nThis is probably due to non-computable boundaries, i.e. no solution for such error allocation.")
  message("Don't need to worry about it.\n")
  return(NaN)
})

return(solution)
}

errs_alpha_only <- function(var_vector, # the variance vector that will be used in calculating information level
                            always_stop_at = NULL, # the stage that the test will always stop, e.g. 3 for H0C
                            alpha = NULL, alphas = NULL, # if want to specify alphas without using error spending approach, set alphas = vector of length K
                            f_err = NULL, # error spending function
                            Imax = NULL, # user-specified Imax
                            ... # ... additional tuning parameters to be passed into f_err
){
  
  if (is.null(alpha) & is.null(alphas)) stop("Either alpha or alphas needs to be specified.")
  
  K <- ifelse(is.null(always_stop_at), length(var_vector), always_stop_at)
  
  I <- rep(NA, K)
  for(k in 1:K){
    I[k] <- 1 / var_vector[k]
  }
  
  if (is.null(Imax)){
    Imax <- I[K]
  }
  
  Iratio <- I/Imax # argument to pass into f_err and g_err
  
  if (is.null(alphas)){
    if (is.null(f_err)) stop("f_err must be specified when alphas is unspecified.")
    fk <- sapply(Iratio, f_err, total_err = alpha, ...)
    fk_ <- sapply(c(0, Iratio[1:(K-1)]), f_err, total_err = alpha, ...)
    alphas <- fk - fk_
    alphas <- c(alphas, rep(0, length(var_vector) - K))
  }
  
  return(alphas)
}

errs <- function(var_vector = NULL, # the variance vector that will be used in calculating information level
                 always_stop_at = NULL, # the stage that the test will always stop, e.g. 3 for H0C
                 alpha = NULL, alphas = NULL, # if want to specify alphas without using error spending approach, set alphas = vector of length K
                 beta = NULL, betas = NULL, # similar as alpha
                 f_err = NULL, g_err = NULL, # error spending functions
                 stg_always_stop_c = NULL, # in computing errors for combined population
                 Imax = NULL, # user-specified Imax
                 ... # ... additional tuning parameters to be passed into f_err and g_err
){
  
  # note: f_err and g_err should take an argument named "total_err"
  
  if (is.null(alpha) & is.null(alphas)) stop("Either alpha or alphas needs to be specified.")
  if (is.null(beta) & is.null(betas)) stop("Either beta or betas needs to be specified.")
  
  # already specify alphas and betas, so simply return them
  if (is.null(alpha) & is.null(beta)) return (rbind(alphas, betas))
  
  # real work
  K <- ifelse(is.null(always_stop_at), length(var_vector), always_stop_at)
  
  if (!is.null(stg_always_stop_c)) {
    K <- stg_always_stop_c
  }
  
  I <- rep(NA, K)
  for(k in 1:K){
    I[k] <- 1 / var_vector[k]
  }
  
  if (is.null(Imax)){
    Imax <- I[K]
  }
  
  Iratio <- I/Imax # argument to pass into f_err and g_err
  
  if (is.null(alphas)){
    if (is.null(f_err)) stop("f_err must be specified when alphas is unspecified.")
    fk <- sapply(Iratio, f_err, total_err = alpha, ...)
    fk_ <- sapply(c(0, Iratio[1:(K-1)]), f_err, total_err = alpha, ...)
    alphas <- fk - fk_
    alphas <- c(alphas, rep(0, length(var_vector) - K))
  }
  if (is.null(betas)){
    if (is.null(g_err)) stop("g_err must be specified when betas is unspecified.")
    gk <- sapply(Iratio, g_err, total_err = beta, ...)
    gk_ <- sapply(c(0, Iratio[1:(K-1)]), g_err, total_err = beta, ...)
    betas <- gk - gk_
    betas <- c(betas, rep(0, length(var_vector) - K))
  }
  
  return(rbind(alphas, betas))
}


# ESF: error spending functions
power_ESF <- function(x, total_err, rho){
  if (x < 0){
    f <- 0
  } else if (x < 1){
    f <- total_err * x^rho
  } else{
    f <- total_err
  }
  return(f)
}

gen_power_ESF <- function(rho){
  function(x, total_err){
    power_ESF(x, total_err, rho)
  }
}

.f_err <- gen_power_ESF(2)
.g_err <- gen_power_ESF(2)