rm(list = ls())

setwd("~/Adaptive-Enrichment-Designs-with-Delayed-Outcomes/for_Josh/analysis")
setwd("../../results_simulation_paper_evaluated/")


generic_filename <- "AD47_pW_wfp_result_pW"
param_name <- "pW"
param_plotlabel <- expression('R'[W]^2)

files_all <- list.files()
files_all <- files_all[grepl(paste0(generic_filename, '*'), files_all)]

result_lists <- list()
info_lists <- list()
for (est in c("ltmle", "unadj")){
  files <- files_all[grepl(paste0('*', est, '*'), files_all)]
  params <- c()
  print(files) # check two sets of params (ltmle, unadj) are exactly the same
  
  temp_result <- list()
  temp_info <- list
  for (case_id in 1:length(files)){
    load(files[case_id])
    params <- c(params, trial_info[[param_name]])
    temp_result <- c(temp_result, list(result))
    temp_info <- c(temp_info, list(trial_info)) 
  }
  print(params)
  
  temp_result <- temp_result[order(params)]
  temp_info <- temp_info[order(params)]
  result_lists[[est]] <- temp_result
  info_lists[[est]] <- temp_info
}
params <- sort(params)

# Compute R-square --------------------------------------------------------

Rsq_W <- readRDS("../dgm/Rsq_empirical.RDS")
params <- Rsq_W[seq(1,length = 10, by = 2), "comb pop_W34"]

# W_to_use <- c("W3", "W4")
# Rsq_W <- sapply(designs_pW_noL, function(l) {
#   tmp <- rep(-9999, 2)
#   for (i in 1:2) {
#     tmp[i] <- t(l$beta.W[W_to_use, i]) %*% .varW[[i]][W_to_use, W_to_use] %*% l$beta.W[W_to_use, i] / .varY["noL", i]
#   }
#   return(tmp)
# })
# 
# Rsq_W <- t(Rsq_W)
# colnames(Rsq_W) <- c("subpop1", "subpop2")
# rownames(Rsq_W) <- paste0("pW", params)
# 
# Rsq_xaxis <- (Rsq_W[, "subpop1"] + Rsq_W[, "subpop2"]) / 2 
# params <- Rsq_xaxis


# check type I error
sapply(result_lists$ltmle, function(r) r$typeIerror)
sapply(result_lists$unadj, function(r) r$typeIerror)

# check power
sapply(result_lists$ltmle, function(r) r$power)
sapply(result_lists$unadj, function(r) r$power)

# check ESS
sapply(result_lists$ltmle, function(r) r$expected_sample_size$ltmle[1, 1:4])
sapply(result_lists$unadj, function(r) r$expected_sample_size$unadj[1, 1:4])

# check duration
sapply(result_lists$ltmle, function(r) r$average_duration$ltmle[1, 1:4])
sapply(result_lists$unadj, function(r) r$average_duration$unadj[1, 1:4])


# Making plots ------------------------------------------------------------

library("reshape2")
library("ggplot2")

setwd("../../plots/AD47")
source("../../src/myggplot.R")



##### plot of power #####

powers <- list()
for (est in c("ltmle", "unadj")){
  powers[[est]] <- sapply(result_lists[[est]], function(r) r$power)
  powers[[est]] <- data.frame(t(powers[[est]]))
  powers[[est]]$params<- params
  powers[[est]]$estimator <- est
}
powers <- rbind(powers$ltmle, powers$unadj)

df_melt <- melt(powers, id = c("params", "estimator"))
ggplot(df_melt, aes(x = params, y = value, col = variable, linetype = estimator)) + 
  geom_line() + 
  xlab(param_plotlabel) +
  ylab("Power") +
  coord_cartesian(xlim = c(0, 0.65), ylim = c(0.75, 1)) +
  scale_colour_manual(values = cbbPalette, guide = FALSE) +
#   scale_color_discrete(guide = FALSE) +
  scale_linetype_discrete(guide = FALSE) +
  #   estimator_legend +
  #   scenario_legend +
  myggfont()

ggsave("pWnoL_power.pdf", height = 4, width = 3.5)



##### plot of ESS #####

ESS <- list()
for (est in c("ltmle", "unadj")){
  ESS[[est]] <- sapply(result_lists[[est]], function(r) r$expected_sample_size[[est]][1, 1:4])
  ESS[[est]] <- data.frame(t(ESS[[est]]))
  ESS[[est]]$params<- params
  ESS[[est]]$estimator <- est
}
ESS <- rbind(ESS$ltmle, ESS$unadj)

df_melt <- melt(ESS, id = c("params", "estimator"))
df_melt$variable <- factor(df_melt$variable, levels = c("H01H02_b", "H11H02_b", "H01H12_b", "H11H12_b"))
ggplot(df_melt, aes(x = params, y = value, col = variable, linetype = estimator)) + 
  geom_line() + 
  xlab(param_plotlabel) +
  ylab("ESS") +
  coord_cartesian(xlim = c(0, 0.65), ylim = c(600, 1550)) +
#   scale_color_discrete(guide = FALSE) +
  scale_colour_manual(values = cbbPalette, guide = FALSE) +
  scale_linetype_discrete(guide = FALSE) +
  #   estimator_legend +
  #   scenario_legend +
  myggfont()
ggsave("pWnoL_ESS.pdf", height = 4, width = 3.5)

# plot of only H01H02 and H11H12
ESS$H11H02_b <- ESS$H01H12_b <- NULL
df_melt <- melt(ESS, id = c("params", "estimator"))
df_melt$variable <- factor(df_melt$variable, levels = c("H01H02_b", "H11H12_b"))
ggplot(df_melt, aes(x = params, y = value, col = variable, linetype = estimator)) + 
  geom_line(size = 2) + 
  xlab(param_plotlabel) +
  ylab("Expected sample size") +
  coord_cartesian(xlim = c(0, 0.65), ylim = c(600, 1550)) +
#   scale_color_discrete(guide = FALSE) +
  scale_colour_manual(values = cbbPalette, guide = FALSE) +
  scale_linetype_discrete(guide = FALSE) +
  #   estimator_legend +
  #   scenario_legend +
  myggfont() +
  blank_theme
ggsave("pWnoL_ESS_JSM.pdf", height = 4, width = 3.5)



##### plot of average duration #####

duration <- list()
for (est in c("ltmle", "unadj")){
  duration[[est]] <- sapply(result_lists[[est]], function(r) r$average_duration[[est]][1, 1:4])
  duration[[est]] <- data.frame(t(duration[[est]]))
  duration[[est]]$params<- params
  duration[[est]]$estimator <- est
}
duration <- rbind(duration$ltmle, duration$unadj)
  
df_melt <- melt(duration, id = c("params", "estimator"))
df_melt$variable <- factor(df_melt$variable, levels = c("H01H02_b", "H11H02_b", "H01H12_b", "H11H12_b"))
ggplot(df_melt, aes(x = params, y = value, col = variable, linetype = estimator)) + 
  geom_line() + 
  xlab(param_plotlabel) +
  ylab("Average duration") +
  coord_cartesian(xlim = c(0, 0.65), ylim = c(4, 9.5)) +
#   scale_color_discrete(guide = FALSE) +
  scale_colour_manual(values = cbbPalette, guide = FALSE) +
  scale_linetype_discrete(guide = FALSE) +
#   estimator_legend +
#   scenario_legend +
  myggfont()
ggsave("pWnoL_duration.pdf", height = 4, width = 3.5)

# plot of only H01H02 and H11H12
duration$H11H02_b <- duration$H01H12_b <- NULL
df_melt <- melt(duration, id = c("params", "estimator"))
df_melt$variable <- factor(df_melt$variable, levels = c("H01H02_b", "H11H12_b"))
ggplot(df_melt, aes(x = params, y = value, col = variable, linetype = estimator)) + 
  geom_line(size = 2) + 
  xlab(param_plotlabel) +
  ylab("Average duration") +
  coord_cartesian(xlim = c(0, 0.65), ylim = c(4, 9.5)) +
#   scale_color_discrete(guide = FALSE) +
  scale_colour_manual(values = cbbPalette, guide = FALSE) +
  scale_linetype_discrete(guide = FALSE) +
  #   estimator_legend +
  #   scenario_legend +
  myggfont() +
  blank_theme
ggsave("pWnoL_duration_JSM.pdf", height = 4, width = 3.5)

