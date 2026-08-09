current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
# devtools::install_github("Eva2703/DensityRegression")
library(DensityRegression)
library(data.table)
source("./densities.R")
source("../../R/mixed_density_smooth.R")
source("./simulation_functions.R")
source("./evaluate_sim_help_functions.R")
library(mgcv)

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

source("evaluate_sim_help_functions.R")

n_simulation_runs <- 200 # number of simulation runs
path <- "./Objects/"
spline_order <- 4
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

scenario_dimensions_mse <- c(sapply(list(scenarios$n_obs, scenarios$step_size), length), 2, 5, 200) # 2 for penalizations, 5 for covariable types, 200 simulation runs
scenario_dimensions_coverage <- c(sapply(list(scenarios$n_obs, scenarios$step_size, covariance_types), length), 2, 5, 200) # 2 for penalizations, 5 for covariable types, 200 simulation runs

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

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
                      n_knots = 10))

covariance_types <- c("Vc", "Vp")
n_obs = c(5000, 10000, 50000, 100000)
n_bins = c(100, 1000, 2000)
penalizations = c("penalized", "unpenalized")
library(latex2exp)
effect_names <- c(TeX("$\\hat{\\f}$"), TeX("$\\hat{\\beta}_0$"),
                  TeX("$\\hat{\\beta}_{West\\_East}$"), TeX("$\\hat{\\beta}_{c\\_age}$"),
                  TeX("$\\hat{\\beta}_{c\\_age, West\\_East}$"),
                  TeX("$\\hat{g}(year)$"), TeX("$\\hat{g}_{West\\_East}(year)$"),
                  TeX("$\\hat{g}_{c\\_age}(year)$"), TeX("$\\hat{g}_{c\\_age, West\\_East}(year)$"))

index_n_obs_step_size <- cbind(rep(seq_along(scenarios$n_obs), each = length(scenarios$step_size)),
                               rep(seq_along(scenarios$step_size), length(scenarios$n_obs)))

# first I reproduce the boxplots of the coverage rates
# every figure contains the effects the number of bins and the number of observations
# this leads to a loop over penalized, Variance type (Vc or Vp) and the different base densities

load_density_scenario <- function(density_params) {
  sim_results <- list()
  coverage_rates <- list()
  for (penalization in penalizations) {
    base_path <- paste0("./simulation_results/", penalization)
    scenario_name <- get_scenario_name(density_params)
    scenario_path <- paste0(base_path, "/", scenario_name)
    save_path <- paste0(scenario_path, "/Simulation_Objects")
    for (N in n_obs) {
      for (G in n_bins) {
        results_path <- paste0(scenario_path, "/n_runs", n_simulation_runs, "_nobs",
                               N, "_nbins", G, ".rds")
        results <- readRDS(results_path)
        sim_results[[penalization]][[paste0("n_obs_", N)]][[paste0("n_bins_", G)]] <- results
      }
      # load the results of the simulations
      coverage_rates[[penalization]]<- readRDS(paste0(save_path, "/coverage_rates.rds"))
    }
  }
  rel_mse_array <- array(numeric(prod(scenario_dimensions_mse)),
                         dim = scenario_dimensions_mse)
  for (l in seq_along(penalizations)) {
    for (i in seq_along(n_obs)) {
      for (j in seq_along(n_bins)) {
        for (n in seq_len(n_simulation_runs)) {
          for (m in seq_len(5)){
            rel_mse_array[i, j, l, m, n] <- sim_results[[penalizations[l]]][[paste0("n_obs_", n_obs[i])]][[paste0("n_bins_", n_bins[j])]][[n]][["MSE_unique_cov"]][["relMSE"]][[m]]
          }
        }
      }
    }
  }
  return (list("coverage_rates" = coverage_rates, "rel_mse_arr" = rel_mse_array, results = sim_results))
}

results_beta_2_2 <- load_density_scenario(list(density_name = "beta", a = 2, b = 2, n_knots = 10))
t <- results_beta_2_2$results$penalized$n_obs_50000
relMSE <- results_beta_2_2$rel_mse_arr
cov_rates <- results_beta_2_2$coverage_rates$penalized
## quick analyisis of coveragee rates
cov_rates[1,3,2,1:5]
# TODO: - create boxplot function for relMSEs
# -  create boxplot function for coverage rate
# - split MSE unique und obs für whole density
# - reorder to match as plot

ymax_main <- 3.8
# To add number of outliers, we cannot plot the whole range
# sort(relMSE[[3]][[1]][4,], decreasing = TRUE)[1:20]
ycut_main <- 3.35
main_effects <- c(1:5)
N_main_effects <- length(main_effects)
params_main <- get_params_for_matrix_plot(n_cols = length(n_bins), n_rows = N_main_effects,
                                          byrow = FALSE, up = 3, le_ri = 8.5)

for (penalization in seq_along(penalizations)) {
  #pdf(paste0("./Images/relMSE_main_", penalizations[penalization], ".pdf"), height = 4, width = 6.5)
  layout(matrix(1:(length(n_bins) * 5), ncol = length(n_bins)),
         heights = c(1, rep(0.47, 5 - 2), 1),
         widths = c(1, rep(0.56, length(n_bins) - 2), 1))
  sapply(seq_along(n_bins),
         function(i)
           sapply(seq_len(N_main_effects),
                  function(k) {
                    param_ind <- (i - 1) * N_main_effects + k
                    if (k == 1) {
                      main <- paste0("G = ", n_bins[i])
                    } else {
                      main <- ""
                    }
                    relMSEs <- lapply(seq_along(n_obs),
                                      function(j) relMSE[j, i, penalization, main_effects[k],])
                    out_of_lim <- lapply(relMSEs, function(r) length(which(r > ycut_main)))
                    relMSEs_plot <- lapply(relMSEs,
                                           function(r) r[which(r <= ycut_main)])
                    par(mar = params_main$mar[[param_ind]])
                    boxplot(relMSEs_plot[length(n_obs):1],
                            # main = main,
                            horizontal = TRUE, lwd = 0.6,
                            # main = paste0("G = ", G[i], ", N = ", N[j]),
                            # ylim = xlims[[i]], # ylim = ylims[[j]],
                            ylim = c(0, ymax_main),
                            xaxt = params_main$xaxt[param_ind], yaxt = "n"
                    )
                    abline(v = ycut_main, lty = 1, lwd = 0.6)
                    sapply(seq_along(out_of_lim), function(n) {
                      if (out_of_lim[[n]] > 0) {
                        text(x = ymax_main * 0.95, y = (length(out_of_lim):1)[n],
                             labels = paste0("+", out_of_lim[n])) # , col = "red")
                      }
                    })
                    mtext(text = main, side = 3, line = 1)
                    if (i == length(n_bins)) {
                      # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                      mtext(side = 4, at = length(n_obs):1, text = paste0("N = ", n_obs),
                            las = 1, line = 1)
                    }
                    if (i == 1) {
                      mtext(text = effect_names[main_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                            side = 2, las = 1, line = 1)
                    }
                    # if (k == N_main_effects) {
                    #   segments(x0 = ycut_main, y0 = 0.5, x1 = ycut_main, y1 = -1.1,
                    #            lwd = 0.6, xpd = TRUE)
                    # }
                  }
           ))
  at_y <- 3 + 3 * (N_main_effects / 2 - 1) * 1.1
  at_x <- - (length(n_bins) / 2 - 1) * ymax_main * 1.1
  mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
  mtext(text = TeX("$relMSE(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
  dev.off()
}
#pdf("./Images/relMSE_main.pdf", height = 4, width = 6.5)
layout(matrix(1:(length(n_bins) * 5), ncol = length(n_bins)),
       heights = c(1, rep(0.47, 5 - 2), 1),
       widths = c(1, rep(0.56, length(n_bins) - 2), 1))
sapply(seq_along(relMSE),
       function(i)
         sapply(seq_len(N_main_effects),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  relMSEs <- lapply(seq_along(N),
                                    function(j) relMSE[[i]][[j]][main_effects[k],])
                  out_of_lim <- lapply(relMSEs, function(r) length(which(r > ycut_main)))
                  relMSEs_plot <- lapply(relMSEs,
                                         function(r) r[which(r <= ycut_main)])
                  par(mar = params_main$mar[[param_ind]])
                  boxplot(relMSEs_plot[length(N):1],
                          # main = main,
                          horizontal = TRUE, lwd = 0.6,
                          # main = paste0("G = ", G[i], ", N = ", N[j]),
                          # ylim = xlims[[i]], # ylim = ylims[[j]],
                          ylim = c(0, ymax_main),
                          xaxt = params_main$xaxt[param_ind], yaxt = "n"
                  )
                  abline(v = ycut_main, lty = 1, lwd = 0.6)
                  sapply(seq_along(out_of_lim), function(n) {
                    if (out_of_lim[[n]] > 0) {
                      text(x = ymax_main * 0.95, y = (length(out_of_lim):1)[n],
                           labels = paste0("+", out_of_lim[n])) # , col = "red")
                    }
                  })
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N),
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[main_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  # if (k == N_main_effects) {
                  #   segments(x0 = ycut_main, y0 = 0.5, x1 = ycut_main, y1 = -1.1,
                  #            lwd = 0.6, xpd = TRUE)
                  # }
                }
         ))
at_y <- 3 + 3 * (N_main_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * ymax_main * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$relMSE(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()


