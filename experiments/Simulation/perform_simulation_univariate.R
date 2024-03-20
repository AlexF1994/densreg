################################################################################
########################## Penalized vs. unpenalized ###########################
################################################################################

# how to structure simulation
# set params -> how to handle different densities -> for now just add in script
# get true density and clr function as well as norm from factory
# prepare simulation
## get approximated density and norm
# run simulation
## sample from simulated density and true density
## fit regression
## calculate coverage probs
## save all estimated densities



current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("./simulation_functions.R")
source("./densities.R")
library(mgcv)

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

scenario_dimensions <- sapply(list(scenarios$n_obs, scenarios$step_size, covariance_types), length)

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

simulate <- function(){

for (penalized in c(FALSE, TRUE)) {
  if (penalized) {
    sp <- NULL
    identifier <- "penalized"
  }
  else {
    sp <- 1
    identifier <- "unpenalized"
  }

  base_path <- paste0("./", identifier)

  for (density_params in scenarios$densities) {
    scenario_name <- get_scenario_name(density_params)
    scenario_path <- paste0(base_path, "/", scenario_name)
    save_path <- paste0(scenario_path, "/Simulation_Objects")

    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

    densities <- get_densities(density_params, calculate_norm = TRUE)
    n_bins <- sapply(scenarios$step_size, function(s) length(seq(0, 1, by = s))) - 1
    theta <- get_theta(densities$clr_density_function)
    saveRDS(theta, paste0(save_path, "/theta.rds"))

    coverage_rate <- list()
    coverage_rate[[scenario_name]] <- array(numeric(prod(scenario_dimensions)),
                                            dim = scenario_dimensions,
                                dimnames = list(paste0("N: ", scenarios$n_obs),
                                                paste0("Bins: ", n_bins),
                                                paste0("Type: ", covariance_types)))

    ##### Perform simulation
    set.seed(scenarios$seed)
    count <- 1
    for(i in seq_along(scenarios$n_obs)) {
      for(j in seq_along(scenarios$step_size)) {
        print(paste(scenario_path, "; G: ", n_bins[j], "; N: ", N_vec[i], "; Count: ", count))
        sim_result <- lapply(1:n_simulation_runs, run_simulation, sample_type = "multinomial",
                             theta = theta, knots = density_params$n_knots, N = scenarios$n_obs[i],
                             step_size = scenarios$step_size[j], alpha = alpha, sp = sp,
                             norm_true = densities$norm_clr_density)
        saveRDS(sim_result, paste0(folder[f], "/Simulation_Objects/S", S, "_N",
                                   N_vec[i], "_G", n_bins[j], ".rds"))
        coverage_rate[[scenario_name]][i, j, 1] <- sum(sapply(1:n_simulation_runs,
                                                              function(l) sim_result[[l]]$coverage_Vc)) / n_simulation_runs
        coverage_rate[[scenario_name]][i, j, 2] <- sum(sapply(1:n_simulation_runs,
                                                              function(l) sim_result[[l]]$coverage_Vp)) / n_simulation_runs
        count <- count + 1
      }
    }
    saveRDS(coverage_rate[[scenario_name]], paste0(folder[f], "/Simulation_Objects/coverage_rates.rds"))
  }
  }
}
################################################################################
########################## Unpenalized with smaller N ##########################
################################################################################

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("./simulation_functions.R")
library(mgcv)

grid <- seq(0, 1, length.out = 100)
a <- c(3, 2, 0.1) #  c(3, 2, 2, 0.1, 3, 2, 0.1)
b <- c(3, 5, NA) # c(3, 5, 5, NA, 3, 5, NA)
K_T <- c(10, 10, 10) # c(10, 10, 30, 10, 10, 10, 10)
ord <- 4
k <- lapply(K_T, function(K) get_knots(K_T = K, ord = ord))
seed <- rep(1542, 3) # c(1542, 1542, 1211, rep(1542, 4))
sp <- list(0, 0, 0) # list(NULL, NULL, NULL, NULL, 0, 0, 0)
scenarios <- c("a3_b3", "a2_b5", "line_a01")
# scenarios <- c("a3_b3", "a2_b5", "a2_b5_K30", "line_a01",
#                "a3_b3", "a2_b5", "line_a01")
folder <- paste0("./unpenalized/", scenarios)
# folder <- c(paste0("./penalized/", scenarios[1:4]), paste0("./unpenalized/", scenarios[5:7]))

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

alpha <- 0.05 # 95% confidence

S <- 200 # simulation runs
N_vec <- c(100, 500, 1000, 2500) # c(5000, 10000, 50000, 100000)
step_size_vec <- c(0.01, 0.001, 0.0005)
type_vec <- c("Vc", "Vp")
lengths_vecs <- sapply(list(N_vec, step_size_vec, type_vec), length)
n_bins <- sapply(step_size_vec, function(s) length(seq(0, 1, by = s))) - 1

coverage_rate <- list()

for (f in seq_along(a)) {
  ##### Get true density
  theta <- readRDS(paste0(folder[f], "/Simulation_Objects/theta.rds"))
  norm_true <- integrate(function(x) interpolate_clr(x, theta = theta, knots = k[[f]])^2, lower = 0, upper = 1)$value

  ##### Perform simulation
  coverage_rate[[f]] <- array(numeric(prod(lengths_vecs)), dim = lengths_vecs,
                              dimnames = list(paste0("N: ", N_vec),
                                              paste0("Bins: ", n_bins),
                                              paste0("Type: ", type_vec)))
  set.seed(1326)
  count <- 1
  for(i in seq_along(N_vec)) {
    for(j in seq_along(step_size_vec)) {
      print(paste(folder[f], "; G: ", n_bins[j], "; N: ", N_vec[i], "; Count: ", count))
      sim_result <- lapply(1:S, run_simulation, sample_type = "multinomial",
                           theta = theta, knots = k[[f]], N = N_vec[i],
                           step_size = step_size_vec[j], alpha = alpha, sp = sp[[f]],
                           norm_true = norm_true)
      saveRDS(sim_result, paste0(folder[f], "/Simulation_Objects/S", S, "_N",
                                 N_vec[i], "_G", n_bins[j], ".rds"))
      coverage_rate[[f]][i, j, 1] <- sum(sapply(1:S, function(l) sim_result[[l]]$coverage_Vc)) / S
      coverage_rate[[f]][i, j, 2] <- sum(sapply(1:S, function(l) sim_result[[l]]$coverage_Vp)) / S
      count <- count + 1
    }
  }
  saveRDS(coverage_rate[[f]], paste0(folder[f], "/Simulation_Objects/coverage_rates_small_N.rds"))
}

################################################################################
######################### Adaptive smoothing parameter #########################
################################################################################

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("./simulation_functions.R")
library(mgcv)

grid <- seq(0, 1, length.out = 100)
a <- c(3, 2, 2, 0.1)
b <- c(3, 5, 5, NA)
K_T <- c(10, 10, 30, 10)
ord <- 4
k <- lapply(K_T, function(K) get_knots(K_T = K, ord = ord))
seed <- c(1542, 1542, 1211, 1542)
sp <- list(NULL, NULL, NULL, NULL)
scenarios <- c("a3_b3", "a2_b5", "a2_b5_K30", "line_a01")
folder <- c(paste0("./penalized/", scenarios[1:4]))

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

alpha <- 0.05 # 95% confidence

S <- 200 # simulation runs
N_vec <- c(5000, 10000, 50000, 100000)
step_size_vec <- c(0.01, 0.001, 0.0005)
type_vec <- c("Vc", "Vp")
lengths_vecs <- sapply(list(N_vec, step_size_vec, type_vec), length)
n_bins <- sapply(step_size_vec, function(s) length(seq(0, 1, by = s))) - 1

coverage_rate <- list()

for (f in seq_along(a)) {
  ##### Get true density
  if (!is.na(b[f])) {
    f_0 <- f_0_beta
  } else {
    f_0 <- f_0_line
  }
  set.seed(seed[f])
  t <- get_interpolation_values(knots = k[[f]], ord = ord)
  theta <- get_theta(a = a[f], b = b[f], K_T = K_T[f], ord = ord, knots = k[[f]], t = t)$theta
  saveRDS(theta, paste0(folder[f], "/Simulation_Objects/adaptive_smoothing/theta.rds"))

  norm_true <- integrate(function(x) interpolate_clr(x, theta = theta, knots = k[[f]])^2, lower = 0, upper = 1)$value

  f_0_clr_interp <- interpolate_clr(x = grid, theta = theta, knots = k[[f]])
  saveRDS(f_0_clr_interp, paste0(folder[f], "/Simulation_Objects/adaptive_smoothing/f_0_clr_interp.rds"))
  par(mfrow = c(1, 2))
  plot(grid, f_0_clr_interp, type = "l", ylab = "clr-density", xlab = "t",
       main = "True clr-density",
       ylim = range(f_0_clr_interp, f_0_clr(grid, a = a[f], b = b[f]), finite = TRUE))
  lines(grid, f_0_clr(grid, a = a[f], b = b[f]), col = "grey")
  points(t, f_0_clr(t, a[f], b[f]), pch = 4)
  legend("bottom", legend = c("True clr-density (interpolation)",
                              paste0("Beta(", a[f], ", ", b[f], ")-clr-density")),
         col = 1:2, lty = 1, bty = "n")
  # # Check:
  # bs_theta <- sapply(seq_along(theta), function(i) bs_L20(k[[f]], grid)[, i] * theta[i])
  # f_0_clr_interp2 <- rowSums(bs_theta) #  c(bs_L20(k[[f]], grid) %*% theta)
  # matplot(grid, bs_L20(k[[f]], grid), type = "l", ylim = range(bs_theta, f_0_clr_interp2))
  # matlines(grid, bs_theta)
  # matlines(grid, f_0_clr_interp2)
  # lines(grid, f_0_clr_interp, col = 2, lty = 2)

  f_0_interp <- interpolate(x = grid, theta = theta, knots = k[[f]])
  saveRDS(f_0_interp, paste0(folder[f], "/Simulation_Objects/adaptive_smoothing/f_0_interp.rds"))
  plot(grid, f_0_interp, type = "l", ylab = "density", xlab = "t",
       main = "True density", ylim = c(0, max(f_0_interp)))
  lines(grid, f_0(grid, a = a[f], b = b[f]), col = "grey")
  legend("bottom", legend = c("True density (interpolation)",
                              paste0("Beta(", a[f], ", ", b[f], ")-density")),
         col = 1:2, lty = 1, bty = "n")

  ##### Perform simulation
  coverage_rate[[f]] <- array(numeric(prod(lengths_vecs)), dim = lengths_vecs,
                              dimnames = list(paste0("N: ", N_vec),
                                              paste0("Bins: ", n_bins),
                                              paste0("Type: ", type_vec)))
  set.seed(1326)
  count <- 1
  for(i in seq_along(N_vec)) {
    for(j in seq_along(step_size_vec)) {
      print(paste(folder[f], "; G: ", n_bins[j], "; N: ", N_vec[i], "; Count: ", count))
      sim_result <- lapply(1:S, run_simulation, sample_type = "multinomial",
                           theta = theta, knots = k[[f]], N = N_vec[i],
                           step_size = step_size_vec[j], alpha = alpha, sp = sp[[f]],
                           norm_true = norm_true, bs = "ad")
      saveRDS(sim_result, paste0(folder[f], "/Simulation_Objects/adaptive_smoothing/S", S, "_N",
                                 N_vec[i], "_G", n_bins[j], ".rds"))
      coverage_rate[[f]][i, j, 1] <- sum(sapply(1:S, function(l) sim_result[[l]]$coverage_Vc)) / S
      coverage_rate[[f]][i, j, 2] <- sum(sapply(1:S, function(l) sim_result[[l]]$coverage_Vp)) / S
      count <- count + 1
    }
  }
  saveRDS(coverage_rate[[f]], paste0(folder[f], "/Simulation_Objects/adaptive_smoothing/coverage_rates.rds"))
}

plot.gam
