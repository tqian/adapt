# change font size, legend position
myggfont <- function(legend_pos = NULL,
                     legend_text_size = 18,
                     legend_title_size = 20,
                     axis_text_size = 16,
                     axis_title_size = 20,
                     plot_title_size = 24) {
  ff <- theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size, face="bold"),
              plot.title = element_text(size = plot_title_size))
  if (is.null(legend_pos)) {
    return(ff)
  } else if (legend_pos == "top-left") {
    return(ff + theme(legend.justification = c(0,1), legend.position = c(0,1)))
  } else if (legend_pos == "top-right") {
    return(ff + theme(legend.justification = c(1,1), legend.position = c(1,1)))
  } else if (legend_pos == "bottom-left") {
    return(ff + theme(legend.justification = c(0,0), legend.position = c(0,0)))
  } else if (legend_pos == "bottom-right") {
    return(ff + theme(legend.justification = c(1,0), legend.position = c(1,0)))
  } else {
    stop("legend_pos needs to be NULL, top-left, top-right, bottom-left, or bottom-right.")
  }
}


# change color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#F0E442", "#0072B2", "#CC79A7")
# from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

myggcolor <- function(palette = "Dark2"){
  return(scale_colour_manual(values = cbbPalette))
}


# make plots to share legend
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, c(lapply(plots, function(x)
      x + theme(legend.position="none")), nrow = 1)),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}



# legend in plotting 4 scenarios in legend (for ESS and duration plot)

scenario_legend_owncolor <- scale_colour_manual(name  = "Scenario",
                                       labels = c(expression(paste("a) ", Delta[1], "=0,", Delta[2], "=0")),
                                                  expression(paste("b) ", Delta[1], "=", delta, "," , Delta[2], "=0")),
                                                  expression(paste("c) ", Delta[1], "=0,", Delta[2], "=", delta)),
                                                  expression(paste("d) ", Delta[1], "=", delta, "," , Delta[2], "=", delta))),
                                       values = cbbPalette)

scenario_legend <- scale_colour_discrete(name  = "Scenario",
                                         labels = c(expression(paste("a) ", Delta[1], "=0,", Delta[2], "=0")),
                                                    expression(paste("b) ", Delta[1], "=", delta, "," , Delta[2], "=0")),
                                                    expression(paste("c) ", Delta[1], "=0,", Delta[2], "=", delta)),
                                                    expression(paste("d) ", Delta[1], "=", delta, "," , Delta[2], "=", delta))))

# legend in plotting 3 hypotheses in legend (for power plot)

hypothesis_legend <- scale_colour_discrete(name  = "To reject",
                                           labels = c(expression(paste('H'['00'])),
                                                      expression(paste('H'['01'])),
                                                      expression(paste('H'['02']))))

# legend for estimators (ltmle and unadj)

estimator_legend <- scale_linetype_discrete(name = "Estimator", labels = c("Adjusted", "Unadjusted"))

estimator_legend_bold <- scale_linetype_manual(name = "Estimator",
                                               labels = c("Adjusted", "Unadjusted"),
                                               values=c("solid", "dotted"))

# power labeller for delayWLY

power_labeller_delayWLY <- function(variable, value){
  return(hypothesis_names[value])
}
hypothesis_names <- list(
  'H0C' = expression(paste("Reject ", 'H'['00'], " under (d)")),
  'H01' = expression(paste("Reject ", 'H'['01'], " under (b)")),
  'H02' = expression(paste("Reject ", 'H'['02'], " under (c)"))
)

# ESS labeller for delayWLY
ESS_labeller_delayWLY <- function(variable, value){
  return(scenario_names[value])
}
scenario_names <- list(
  'H01H02_b' = expression(paste("a) ", Delta[1], "=0,", Delta[2], "=0")),
  'H11H02_b' = expression(paste("b) ", Delta[1], "=", delta, "," , Delta[2], "=0")),
  'H01H12_b' = expression(paste("c) ", Delta[1], "=0,", Delta[2], "=", delta)),
  'H11H12_b' = expression(paste("d) ", Delta[1], "=", delta, "," , Delta[2], "=", delta))
)
blank_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

myggsave <- function(filename = default_name(plot), height= 4, width= 4, dpi= 72, ...) {
  ggsave(filename=filename, height=height, width=width, dpi=dpi, ...)
}
