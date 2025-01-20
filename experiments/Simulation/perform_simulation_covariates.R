current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("./simulation_functions.R")
source("./densities.R")
source("../../R/mixed_density_smooth.R")
library(mgcv)

set.seed(25)
alpha <- 0.05 # 95% confidence
n_simulation_runs <- 200 # number of simulation runs
covariance_types <- c("Vc", "Vp")

spline_order <- 4
#sp <- list(NULL, NULL, NULL, NULL, 0, 0, 0) # smoothing parameter (NULL is default, 0 corresponds to no penalization)
scenarios <- list(
  densities = list(list(density_name = "beta", a = 2, b = 2, n_knots = 10),
                   list(density_name = "beta", a = 2, b = 3, n_knots = 10),
                   list(density_name = "beta", a = 3, b = 3, n_knots = 10),
                   list(density_name = "beta", a = 5, b = 5, n_knots = 10),
                   list(density_name = "truncated_normal", mean = 0, sd = 1,
                        n_knots = 10),
                   list(density_name = "truncated_normal", mean = 1, sd = 1,
                        n_knots = 10),
                   list(density_name = "truncated_normal", mean = 1, sd = 2,
                        n_knots = 10),
                   list(density_name = "truncated_normal", mean = 2, sd = 2,
                        n_knots = 10)),
  seed = 1542,
  n_obs = c(5000, 10000, 50000, 100000),
  step_size = c(0.01, 0.001, 0.0005)
)

scenario_dimensions <- c(sapply(list(scenarios$n_obs, scenarios$step_size, covariance_types), length), 5) #  5 for covariable types

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

simulate_with_covariates <- function(){

  for (penalized in c(FALSE, TRUE)) {
    set.seed(scenarios$seed)
    if (penalized) {
      sp <- NULL
      identifier <- "penalized"
    }
    else {
      sp <- 0
      identifier <- "unpenalized"
    }

    base_path <- paste0("./covariates/", identifier)

    for (density_params in scenarios$densities) {
      print(density_params$density_name)
      scenario_name <- get_scenario_name(density_params)
      scenario_path <- paste0(base_path, "/", scenario_name)
      save_path <- paste0(scenario_path, "/Simulation_Objects")

      if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

      # Problem wir samplen immer noch aus der spline approx. der density --> to do
      n_bins <- sapply(scenarios$step_size, function(s) length(seq(0, 1, by = s))) - 1
      n_obs_approx <- 1000
      approx_design_matrix <- sample_covariates(n_obs_approx)
      approx_results <- get_approx_results_with_covariates(density_params = density_params,
                                                           covariates =  approx_design_matrix)
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

      ##### Perform simulation
      count <- 1
      for(i in seq_along(scenarios$n_obs)) {
        X <- sample_covariates(scenarios$n_obs[i])
        for (j in seq_along(scenarios$step_size)) {
          print(paste(scenario_path, "; G: ", n_bins[j], "; N: ", scenarios$n_obs[i], "; Count: ", count))
          sim_result <- lapply(1:n_simulation_runs, run_simulation_with_covariates, sample_type = "bin",
                               approx_results = approx_results, n_knots = density_params$n_knots, n_obs = scenarios$n_obs[i],
                               step_size = scenarios$step_size[j], alpha = alpha, sp = sp, bs = "md", covariates = X)
          saveRDS(sim_result, paste0(scenario_path, "/n_runs", n_simulation_runs, "_nobs",
                                     scenarios$n_obs[i], "_nbins", n_bins[j], ".rds"))
          coverage_rate[[scenario_name]][i, j, 1] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_Vc)) / n_simulation_runs
          coverage_rate[[scenario_name]][i, j, 2] <- sum(sapply(1:n_simulation_runs,
                                                                function(l) sim_result[[l]]$coverage_Vp)) / n_simulation_runs
          count <- count + 1
        }
      }
      saveRDS(coverage_rate[[scenario_name]], paste0(save_path, "/coverage_rates.rds"))
    }
  }
}

