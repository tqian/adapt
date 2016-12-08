### attach enrollment time after simulating a data set ###

## Functions: ##

# attach_time.wrapper: wrapper of attach_time to be used in simu_eval_trial
# attach_time: attach enrollment time to data set, using constant or random arrival time

# generate_Poisson_arrival: generate Poisson arrival time


attach_time.wrapper <- function(dt, setup) {
  
  if (!is.null(setup$enrollmethod)) {
    result <- attach_time(dt,
                          cumss_sub1 = setup$cumss_sub1,
                          cumss_sub2 = setup$cumss_sub2,
                          erate_sub1 = setup$erate_sub1,
                          delay_WL = setup$delay_WL,
                          delay_LY = setup$delay_LY,
                          method = setup$enrollmethod)
  } else {
    result <- attach_time(dt,
                          cumss_sub1 = setup$cumss_sub1,
                          cumss_sub2 = setup$cumss_sub2,
                          erate_sub1 = setup$erate_sub1,
                          delay_WL = setup$delay_WL,
                          delay_LY = setup$delay_LY,
                          method = "Poisson-homo")
  }  
  return(result$dt)  
}


attach_time <- function(dt, cumss_sub1, cumss_sub2, erate_sub1, delay_WL, delay_LY,
                        method = c("continuous", "Poisson-homo")) {
  
  # dt is the data frame provided by Liz, with columns:
  # subpop (subpopulation), W (baseline var), A (trt assignment), L (short-term outcome), Y (final outcome)

  erate_sub2 <- erate_sub1 * cumss_sub2[1] / cumss_sub1[1]
  nmax_sub1 <- cumss_sub1[length(cumss_sub1)]
  nmax_sub2 <- cumss_sub2[length(cumss_sub2)]
  
  # Assign enrollment time, short-term time and final time
  dt$T_enroll <- rep(NA, nrow(dt))
  
  method <- match.arg(method)

  if (method == "continuous") {
    dt$T_enroll[dt$subpop == 1] <- sample.int(nmax_sub1) / erate_sub1
    dt$T_enroll[dt$subpop == 2] <- sample.int(nmax_sub2) / erate_sub2
  } else if (method == "Poisson-homo") {
    dt$T_enroll[dt$subpop == 1] <- generate_Poisson_arrival(nmax_sub1, erate_sub1)
    dt$T_enroll[dt$subpop == 2] <- generate_Poisson_arrival(nmax_sub2, erate_sub2)
  }

  dt$T_shortterm <- dt$T_enroll + delay_WL
  dt$T_final <- dt$T_shortterm + delay_LY
  
  # Calculate interim time according to the cumulative sample (with final outcome available) for subpop 1
  
  interim_time <- sort(dt$T_final[dt$subpop == 1])[cumss_sub1]
  
#   interim_time <- cumss_sub1 / erate_sub1 + delay_WL + delay_LY + .Machine$double.eps^(1/2)
#   # Check interim_time to ensure exact cumulative sample size
#   cumss_sub1_check <- sapply(interim_time, function(time){sum(subset(dt, subpop == 1)$T_final <= time)})
#   cumss_sub2_check <- sapply(interim_time, function(time){sum(subset(dt, subpop == 2)$T_final <= time)})
#   
#   if (!(all(cumss_sub1_check == cumss_sub1) & all(cumss_sub2_check == cumss_sub2))) {
#     stop("Cumulative sample size error! Check interim time and cumss in dgm_cumss.")
#   }
  
  return(list(dt = dt, interim_time = interim_time))
  
}


generate_Poisson_arrival <- function(n, rate) {
  cumsum(rexp(n, rate))
}



