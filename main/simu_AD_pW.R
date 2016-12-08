rm(list = ls())

library(ltmle)
library(mvtnorm)


##### Read in design #####

# first, set working directory to current file location
setwd("~/Adaptive-Enrichment-Designs-with-Delayed-Outcomes/for_Josh/main")
setwd("..")
load("data/ADTianchen_new.RData")

designs_pW <- readRDS("dgm/AD_designs_pW.RDS")

DESIGN_AD <- 47

token_para_result <- "AD47_pW_wfp"
result_gathered <- "gathered_AD47_pW_wfp"
eval_result_filename <- "AD47_pW_wfp_result"

alpha_c <- NULL
alpha_1 <- NULL
alpha_2 <- NULL
alphas_user <- c(0.0027569537, 0.0005642895, 0.0008482876, 0.0012750645, 0.0011654439,
                 0.0006570932, 0.0007052451, 0.0028025699, 0.0014532974, 0.0037450390,
                 0.000112055, 0.002327773, 0.001161325, 0.002593472, 0.002668550)
l_1k <- c(-4.1222874, 0.4032743, -1.4766927, 0.9381067, -3.4339860)
l_2k <- c(-0.09831440, 0.28603462, 0.42349497, 0.92976045, -0.04665651)

dgm_srcfile <- "dgm/AD_dgm_new.R"
dgm_name <- "dgm"

num_simu <- 1000
npara <- 20

betas_user <- NULL
f_err <- NULL

delay_LY <- 1
delay_WL <- 1

W_to_use <- c("W3", "W4")
L_to_use <- 'L'

exL <- TRUE

nmax1 <- 955.4931 * 2
nmax2 <- 1081.69 * 2

###########################
##### Simulation part #####
###########################

# load data
source(dgm_srcfile)
setwd("src")
source("simulation.analysis_fixed_time.R")
setwd("..")

# parallel job id and random seeds
version_fut <- as.integer(Sys.getenv("SGE_TASK_ID"))
if (is.na(version_fut)) {
    version_fut <- npara
}
taskID <- (version_fut %% npara) + 1
parallel_seeds <- data.matrix(read.csv("misc/parallel_seeds.csv"))
.Random.seed <- parallel_seeds[taskID, ]
case_id <- (version_fut - 1) %/% npara + 1

# decide the parameters of current job
current_design <- designs_pW[[case_id]]
pW <- current_design$pW[1]
betas <- current_design$betas
alphas <- current_design$alphas
sigmaY <- current_design$sigmaY
sigmaL <- current_design$sigmaL

# update file names to be saved
appen <- paste0("pW", pW)
token_para_result <- paste0(token_para_result, "_", appen)
result_gathered <- paste0(result_gathered, "_", appen)
eval_result_filename <- paste0(eval_result_filename, "_", appen)

##### import timer and nmaxs
timer_design <- readRDS("design_timer/timers_AD47_pW_wfp.rds")

for (est in c("unadj", "ltmle")) {
  
  print("###############")
  print(est)
  print("###############")
  
  p1 <- nmax1 / (nmax1 + nmax2)
  timer <- timer_design[[est]]$timers[as.character(pW), ]
  
  erate_sub1 <- (500 / 3) * p1
  
  # print out simulation specification
  print(version_fut)
  print(taskID)
  print(nmax1)
  print(nmax2)
  print(timer)
  print(.Random.seed)
  print(case_id)
  print(appen)
  print(pW)
  print(delay_WL)
  print(delay_LY)
  print(token_para_result)
  print(result_gathered)
  print(eval_result_filename)
  
  # record trial_info
  trial_info <- list(num_simu = num_simu, npara = npara,
                     nmax1 = nmax1, nmax2 = nmax2,
                     timer = timer,
                     erate_sub1 = erate_sub1,
                     pW = pW,
                     delay_WL = delay_WL, delay_LY = delay_LY,
                     W_to_use = W_to_use, L_to_use = L_to_use,
                     dgm_srcfile = dgm_srcfile, dgm_name = dgm_name)
  
  nsim <- num_simu / npara
  
  eval(parse(text = paste0("dgm <- ", dgm_name)))
  
  filenames_H01H02 <- paste0("results_simulation_paper/", token_para_result, "_", est, "_H01H02_", 1:npara, ".rda")
  filenames_H11H02 <- paste0("results_simulation_paper/", token_para_result, "_", est, "_H11H02_", 1:npara, ".rda")
  filenames_H01H12 <- paste0("results_simulation_paper/", token_para_result, "_", est, "_H01H12_", 1:npara, ".rda")
  filenames_H11H12 <- paste0("results_simulation_paper/", token_para_result, "_", est, "_H11H12_", 1:npara, ".rda")
  
  # if the Result directory doesn't exist, create it
  dir.create("results_simulation_paper/", showWarnings = FALSE)
  
  for (scenario in c("H01H02", "H11H02", "H01H12", "H11H12")){
    test <- simu_trial_fix_time(nsim = nsim, W_to_use = W_to_use, L_to_use = L_to_use,
                                  dgm = dgm, print_sim_progress = TRUE,
                                  timer = timer, last_timer_stop_enrollment = TRUE,
                                  nmax1 = nmax1, nmax2 = nmax2,
                                  erate_sub1 = erate_sub1, delay_WL = delay_WL, delay_LY = delay_LY,
                                  betas = betas,
                                  alphas = alphas,
                                  sigmaY = sigmaY,
                                  sigmaL  = sigmaL,
                                  data = df, S = "subpop", W = c("W1", "W2", "W3", "W4", "W5"),
                                  AD_dgm_scn = scenario, estimator = est,
                                  enroll_method = "Poisson-homo")
    
    eval(parse(text = paste0("filename = filenames_", scenario, "[taskID]")))
    save(test, file = filename)
  }
  
  
  
  #######################################
  ##### Gather results and Evaluate #####
  #######################################
  
  if (taskID == 1) { # the 1st parallel job will wait for all to complete, then do the gathering
    setwd("src/")
    source("gather_results.R")
    source("evaluation.errsp.analysis_fixed_time.H00H01H02.R")
    setwd("..")
    
    while (!all(all(file.exists(filenames_H01H02)), all(file.exists(filenames_H11H02)),
                all(file.exists(filenames_H01H12)), all(file.exists(filenames_H11H12)))) {
      s0102 <- sum(file.exists(filenames_H01H02))
      s1102 <- sum(file.exists(filenames_H11H02))
      s0112 <- sum(file.exists(filenames_H01H12))
      s1112 <- sum(file.exists(filenames_H11H12))
      cat(paste0("\nWaiting for parallel jobs to be done:\n"))
      cat(paste0("    H01H02: ", s0102, "/", npara, "\n"))
      cat(paste0("    H11H02: ", s1102, "/", npara, "\n"))
      cat(paste0("    H01H12: ", s0112, "/", npara, "\n"))
      cat(paste0("    H11H12: ", s1112, "/", npara, "\n"))
      
      Sys.sleep(30)
    }
    
    Sys.sleep(30)
    
    trials_H01H02 <- gather_results_timer(filenames_H01H02)
    trials_H11H12 <- gather_results_timer(filenames_H11H12)
    trials_H11H02 <- gather_results_timer(filenames_H11H02)  
    trials_H01H12 <- gather_results_timer(filenames_H01H12)
    
    dir.create("results_simulation_paper_gathered/", showWarnings = FALSE)
    save(trial_info, trials_H01H02, trials_H11H12, trials_H11H02, trials_H01H12,
         file = paste0("results_simulation_paper_gathered/", result_gathered, "_", est, ".rda"))
    
    f_err <- gen_power_ESF(2)
    
    result <- eval_trial_fix_time_H00H01H02(trials_H01H02, trials_H11H12, trials_H11H02, trials_H01H12,
                                             nmax1 = nmax1, nmax2 = nmax2, timer = timer,
                                             erate_sub1 = erate_sub1, delay_WL = delay_WL, delay_LY = delay_LY,
                                             alphas_user = alphas_user,
                                             alpha_c = alpha_c, alpha_1 = alpha_1, alpha_2 = alpha_2,
                                             f_err = f_err,
                                             binding_fut = FALSE,
                                             modify_l = "specify-vector", l_1k = l_1k, l_2k = l_2k,
                                             estimator = est, print = "more")
    
#     print(result)
    
    dir.create("results_simulation_paper_evaluated/", showWarnings = FALSE)
    save(result, trial_info, file = paste0("results_simulation_paper_evaluated/",
                                           eval_result_filename, "_", est, ".rda"))
    
    # remove parallel results
    file.remove(filenames_H01H02)
    file.remove(filenames_H11H02)
    file.remove(filenames_H01H12)
    file.remove(filenames_H11H12)
  } # end if (taskID == 1)
  
} # end for loop c("unadj", "ltmle")