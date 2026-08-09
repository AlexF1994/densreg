current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("./densities.R")
source("../../R/mixed_density_smooth.R")
source("./simulation_functions.R")
library(mgcv)
library(parallel)

set.seed(25)
alpha <- 0.05 # 95% confidence
n_simulation_runs <- 200 # number of simulation runs
covariance_types <- c("Vc", "Vp")

spline_order <- 4
scenarios <- list(
  densities = list(list(density_name = "truncated_normal", mean = 0.1, sd = 1,
                        n_knots = 10),
                   list(density_name = "truncated_normal", mean = 0.1, sd = 1,
                        n_knots = 10),
                   list(density_name = "truncated_normal", mean = 0.5, sd = 2,
                        n_knots = 10),
                   list(density_name = "truncated_normal", mean = 0.5, sd = 2,
                        n_knots = 10)),
  seed = 1542,
  n_obs = c(50000, 150000, 500000, 1000000),
  step_size = c(0.02, 0.01, 0.005)
)

scenario_dimensions <- c(sapply(list(scenarios$n_obs, scenarios$step_size, covariance_types), length), 5) #  5 for covariable types

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

index_n_obs_step_size <- cbind(rep(seq_along(scenarios$n_obs), each = length(scenarios$step_size)),
                               rep(seq_along(scenarios$step_size), length(scenarios$n_obs)))

simulate_with_covariates <- function(n_scenario, which, parallel = FALSE) {
  if (parallel) {
    stop("Currently, only parallel execution is supported. Please set 'parallel = TRUE'.")
  # STEP 1: SET UP PARALLEL CLUSTER
  # Detect available cores and leave one free for system processes
  n_cores <- max(1, detectCores() - 1)
  cat("Using", n_cores, "cores for parallel processing\n")

  # Create cluster object
  cl <- makeCluster(n_cores)

  # STEP 2: PREPARE CLUSTER WITH DEPENDENCIES
  # Export necessary functions and variables to all worker processes
  clusterExport(cl, c(
    # Variables
    "alpha", "scenarios",
    # Functions from your source files (add any others you need)
    "run_simulation_with_covariates", "sample_covariates",
    "get_knots", "get_approx_results_with_covariates",
    "sum_constrained_spline_design_matrix", "get_scenario_name"
  ), envir = environment())

  # Load required libraries on each worker
  clusterEvalQ(cl, {
    library(data.table)
    library(mgcv)
    library("r2r")
    library("purrr")
    library("truncnorm")
    library("rlist")
  })

  # Source required files on each worker
  clusterEvalQ(cl, {
    source("./densities.R")
    source("../../R/mixed_density_smooth.R")
    source("./simulation_functions.R")
  })
  } else {
    # If not running in parallel, just print a message
    cat("Running simulations sequentially (not in parallel)\n")
  }

  for (penalized in c(FALSE, TRUE)) {
    set.seed(scenarios$seed)
    if (penalized) {
      sp <- -1
      identifier <- "penalized"
    }
    else {
      sp <- 0
      identifier <- "unpenalized"
    }

    base_path <- paste0("./covariates/", identifier)

    for (n_s in n_scenario) {
      set.seed(scenarios$seed)
      density_params <- scenarios$densities[[n_s]]
      print(density_params$density_name)
      scenario_name <- get_scenario_name(density_params)
      scenario_path <- paste0(base_path, "/", scenario_name)
      save_path <- paste0(scenario_path, "/Simulation_Objects")

      if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

      n_bins <- sapply(scenarios$step_size, function(s) length(seq(0, 1, by = s))) - 1
      n_obs_approx <- 2000
      range_smooth_covariates <- c(1, 5) # avoid 0 as it leads to multicollinearity with base component
      range_linear_covariates <- c(1, 5) # avoid zero as otherwise constant density is introduced
      approx_design_matrix <- sample_covariates(n_obs_approx, range_smooth_covariates, range_linear_covariates)
      # note that you have to change the range if you change the method for sampling smooth covariates
      grid_hist <- seq(from = 0, to = 1, by = 0.01)
      knots_smooth_covariate <- get_knots(n_splines = 8, ord = 4, range_ = range_smooth_covariates)
      quantiles_density <- grid_hist[1:(length(grid_hist) - 1)] + 0.01 / 2
      approx_results <- get_approx_results_with_covariates(density_params = density_params,
                                                           covariates =  approx_design_matrix,
                                                           knots_smooth_covariate = knots_smooth_covariate,
                                                           sp = sp, quantiles_density = quantiles_density)
      saveRDS(approx_results$theta, paste0(save_path, "/theta.rds"))

      coverage_rate <- list()
      coverage_rate[[scenario_name]] <- array(numeric(prod(scenario_dimensions)),
                                              dim = scenario_dimensions,
                                              dimnames = list(paste0("N: ", scenarios$n_obs),
                                                              paste0("Bins: ", n_bins),
                                                              paste0("Type: ", covariance_types),
                                                              paste0("density_component: ",  c("base_density",
                                                                                              "binary_component",
                                                                                              "linear_component",
                                                                                              "smooth_component",
                                                                                              "complete_density"))))
      if (parallel) {
        # STEP 3: EXPORT SCENARIO-SPECIFIC VARIABLES TO CLUSTER
        clusterExport(cl, c("sp", "approx_results", "density_params"), envir = environment())
      }
      # STEP 3: EXPORT SCENARIO-SPECIFIC VARIABLES TO CLUSTER
      ##### Perform simulation
      count <- 1
      for (k in which) {
        set.seed(scenarios$seed)
        i <- index_n_obs_step_size[k, 1]
        X <- sample_covariates(scenarios$n_obs[i], range_smooth_covariates, range_linear_covariates)
        smooth_covariate_design_matrix <- sum_constrained_spline_design_matrix(X$smooth_variable, knots_smooth_covariate)
        if (parallel) {
          # STEP 4: EXPORT LOOP-SPECIFIC VARIABLES TO CLUSTER
          clusterExport(cl, c("X", "smooth_covariate_design_matrix"), envir = environment())
        }
        #for (j in seq_along(scenarios$step_size)) {
        j <- index_n_obs_step_size[k, 2]
          if (parallel) {
            # STEP 5: REPLACE lapply() WITH parLapply() FOR PARALLEL EXECUTION
            # This is the main parallelization - instead of running 200 simulations sequentially,
            # we distribute them across multiple CPU cores
            sim_result <- parLapply(cl, 1:n_simulation_runs, function(run_index) {
              # Set different seed for each simulation run to ensure reproducibility
              # but still have variation between runs
              set.seed(scenarios$seed + run_index)

              run_simulation_with_covariates(
                run_index,
                sample_type = "bin",
                approx_results = approx_results,
                n_knots = density_params$n_knots,
                n_obs = scenarios$n_obs[i],
                step_size = scenarios$step_size[j],
                alpha = alpha,
                sp = sp,
                bs = "md",
                covariates = X,
                smooth_covariate_design_matrix = smooth_covariate_design_matrix
              )
            })
          } else {
            print(paste(scenario_path, "; G: ", n_bins[j], "; N: ", scenarios$n_obs[i], "; Count: ", count))
            # If not running in parallel, use regular lapply
            sim_result <- lapply(1:n_simulation_runs, run_simulation_with_covariates, sample_type = "bin",
                                 approx_results = approx_results, n_knots = density_params$n_knots, n_obs = scenarios$n_obs[i],
                                 step_size = scenarios$step_size[j], alpha = alpha, sp = sp, bs = "md", covariates = X,
                                 smooth_covariate_design_matrix = smooth_covariate_design_matrix)
          }
          saveRDS(sim_result, paste0(scenario_path, "/n_runs", n_simulation_runs, "_nobs",
                                     scenarios$n_obs[i], "_nbins", n_bins[j], ".rds"))

          # base component
          coverage_rate[[scenario_name]][i, j, 1, 1] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_base$check_coverage_Vc)) / n_simulation_runs
          coverage_rate[[scenario_name]][i, j, 2, 1] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_base$check_coverage_Vp)) / n_simulation_runs
          # binary component
          coverage_rate[[scenario_name]][i, j, 1, 2] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_binary$check_coverage_Vc)) / n_simulation_runs
          coverage_rate[[scenario_name]][i, j, 2, 2] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_binary$check_coverage_Vp)) / n_simulation_runs

          # linear component
          coverage_rate[[scenario_name]][i, j, 1, 3] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_linear$check_coverage_Vc)) / n_simulation_runs
          coverage_rate[[scenario_name]][i, j, 2, 3] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_linear$check_coverage_Vp)) / n_simulation_runs

          # smooth component
          coverage_rate[[scenario_name]][i, j, 1, 4] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_smooth$check_coverage_Vc)) / n_simulation_runs
          coverage_rate[[scenario_name]][i, j, 2, 4] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_smooth$check_coverage_Vp)) / n_simulation_runs

          # whole density
          coverage_rate[[scenario_name]][i, j, 1, 5] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_whole_density$check_coverage_Vc)) / n_simulation_runs
          coverage_rate[[scenario_name]][i, j, 2, 5] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_whole_density$check_coverage_Vp)) / n_simulation_runs
          count <- count + 1
        #}
      }
      saveRDS(coverage_rate[[scenario_name]], paste0(save_path, "/coverage_rates.rds"))
    }
  }
  if (parallel) {
    # STEP 6: CLEAN UP - ALWAYS STOP THE CLUSTER WHEN DONE
    stopCluster(cl)
    cat("Parallel cluster stopped\n")
  }
}

