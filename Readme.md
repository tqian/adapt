# Simulate Adaptive Enrichment Design

This code simulates adaptive enrichment design with 2 subpopulations. The multiple testing procedures used is detailed in: Rosenblum M., Qian T., Du Y., and Qiu H. (2016) Multiple Testing Procedures for Adaptive Enrichment Designs: Combining Group Sequential and Reallocation Approaches. Biostatistics, 17(4), 650–662.




## Code Structure

main/: code that can be run in parallel.

    estvar: estimate variance at a fine time grid (to determine analysis timing)

    simu: simulate trials and evaluate trials (i.e. compute type I error, power, etc)

src/: source code of functions that are used in main.

    simulation.analysis_fixed_time.R: simulate trials, when each interim analysis takes place at pre-determined time.

    evaluation.errsp.analysis_fixed_time.H00H01H02.R: evalute trial (determine the outcome of each trial, compute type I error, power, etc) using error spending approach, when each interim analysis takes place at pre-determined time.

    gather_results.R: gather parallelly simulated trials (so that estimators form a huge matrix), and computes the covariance matrix of estimators at different stages / for different subpopulations.

dgm/: source code of data generating mechanism

data/: folder containing the data set

design_timer/:
    estimate_variance_wfp.R: code that uses the result from main/estvar to determine the time of each interim anslysis (i.e. conduct interim analysis once variance hit certain threshold). wfp means "wait-for-pipeline".

analysis/: folder containing code to make plots

## Folders that will be generated by running code in main/:

results_simulation_paper/: temporary results from parallel simulation.

results_simulation_paper_estvar/: mean and covariance matrix of estimators for "estvar"

results_simulation_paper_gathered/: trial estimators for "simu" and "estvar" (gathered from parallel simulations)

results_simulation_paper_evaluated/: evaluated trial results, containing type I error, power, etc.