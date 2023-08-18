################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size
# options(scipen = scipen)

path <- "./Simulation_Images/"

alpha <- 0.05
S <- 200 # simulation runs
N_vec <- c(5000, 10000, 50000, 100000)
step_size_vec <- c(0.01, 0.001, 0.0005)
n_bins <- sapply(step_size_vec, function(s) length(seq(0, 1, by = s))) - 1

N_long <- rep(N_vec, each = length(n_bins))
n_bins_long <- rep(n_bins, length(N_vec))
sim_result <- lapply(1:length(N_long), function(i)
  readRDS(paste0("./Simulation_Objects/S", S, "_sample_size", N_long[i], "_bins", n_bins_long[i], ".rds")))
names(sim_result) <- paste0("Sample Size: ", N_long, ", Bins: ", n_bins_long, ", Simulation runs: ", S)

coverage_rates <- readRDS("./Simulation_Objects/coverage_rates.rds")
f_0_interp <- readRDS("./Simulation_Objects/f_0_interp.rds")
f_0_clr_interp <- readRDS("./Simulation_Objects/f_0_clr_interp.rds")
grid <- seq(0, 1, length.out = nrow(f_0_interp))

################################################################################
############################# Plot coverage rates ##############################
################################################################################

# Coverage rates according to our approach
library(plot.matrix)
plot_results <- function(result, S = 200, type = c("Vc", "Vp")) {
  plot(result, xlab = "", ylab = "Sample size", axis.col = NULL,
       axis.row = NULL, main = paste0("Coverage rates in ", S, " simulation runs using ", type))
  axis(1, at = 1:3, labels = n_bins)
  axis(2, at = 4:1, labels = N_vec, las = 2)
  text(rep(1:3, each = 4), rep(1:4, 3), labels = result[nrow(result):1, ])
  mtext("Bin number", side = 1, line = 3, at = dim(result)[2] / 2 + 0.5)
}

def.par <- par(no.readonly = TRUE)
pdf(paste0(path, "Simulation_CI_200.pdf"), height = 9)
par(mar = c(mar = c(5, 4.7, 4, 4) + 0.1), mgp = c(3.7, 1, 0), mfrow = c(2, 1))
plot_results(coverage_rates[, , 1], S = 200, type = "Vc")
plot_results(coverage_rates[, , 2], S = 200, type = "Vp")
dev.off()
# par(def.par)

# Coverage rage averaged across the domain of the function (Wood, 2017, Section 6.10/p. 294)
coverage_pw <- lapply(seq_along(N_long),
                      function(i) sapply(1:S,
                                         function(j) sum(sim_result[[i]][[j]]$coverage_pw$CI_p_check) / n_bins_long[i]))

names(coverage_pw) <- paste0("Sample Size: ", N_long, ", Bins: ", n_bins_long) #, ", Simulation runs: ", S_long)

pdf(paste0(path, "Simulation_CI_average_across_function.pdf"), height = 9, width = 8)
par(mfrow = c(4, 3), mar = c(1, 4, 4, 1) + 0.1)
sapply(1:length(N_long), function(i) {
  boxplot(coverage_pw[[i]], main = names(coverage_pw)[[i]], ylim = c(0, 1))
  abline(h = 1 - alpha, col = "green")
  abline(h = mean(coverage_pw[[i]]), col = "red")
  })
dev.off()
# par(def.par)

################################################################################
###################### Plot fit vs. histograms and truth #######################
################################################################################
# For each combination of N (sample size of individual observations) and n_bins
# (number of histogram bins): Sort results by test statistic value
statistics <- lapply(seq_along(N_long),
                     function(n) {
                       r <- cbind(Vc = sapply(1:S, function(s) sim_result[[n]][[s]]$statistic_vc),
                                  Vp = sapply(1:S, function(s) sim_result[[n]][[s]]$statistic_vp))
                       rownames(r) <- 1:S
                       r
                     })
sim_result_sort_Vc <- lapply(seq_along(N_long), function(n) sim_result[[n]][order(statistics[[n]][, 1])])
sim_result_sort_Vp <- lapply(seq_along(N_long), function(n) sim_result[[n]][order(statistics[[n]][, 2])])
names(sim_result_sort_Vc) <- names(sim_result_sort_Vp) <- names(sim_result)

sim_result_dat <- sim_result_sort_Vc # Can also be replaced with sorting by Vp or no sorting
# I tried to add a drop-down menu for changing the data set (sorting), but the
# resulting manipulate-plot is super slow... Thus, the user has to chose the data
# set to use (sim_result_dat) manually before running manipulate.
library(manipulate)
manipulate({
  par(mfrow = c(1, 2))
  # # Using picker "type" to change data sets works, but is super slow
  # if (type == "Sort by Vc") {
  #   sim_result_dat <- sim_result_sort_Vc
  # } else if (type == "Sort by Vp") {
  #   sim_result_dat <- sim_result_sort_Vp
  # } else {
  #   sim_result_dat <- sim_result
  # }
  n <- which(N_long == N & n_bins_long == n_bins)
  dta_dens <- sim_result_dat[[n]][[s]]$data
  step_size <- max(diff(dta_dens$y))
  main <- paste0("Sample Size: ", N_long[n], ", Bins: ", n_bins_long[n],
                 ", Iteration number: ", order(statistics[[n]][, 1])[s])
  # Bayes level
  hist(rep(dta_dens$y, dta_dens$counts), breaks = seq(0, 1, by = step_size),
       right = FALSE, xlab = "y", ylab = "density", freq = FALSE, main = "")
  lines(grid, f_0_interp, col = "green")
  lines(dta_dens$y, sim_result_dat[[n]][[s]]$f_hat, col = "red")
  # clr-level
  plot(grid, f_0_clr_interp, col = "green", type = "l", xlab = "y", ylab = "clr-density",
       ylim = range(sapply(1:S, function(l) range(sim_result_dat[[n]][[l]]$f_hat_clr)),
                    f_0_clr_interp))
  mtext(main, line = 1, at = -0.4, cex = 1.3)
  lines(dta_dens$y, sim_result_dat[[n]][[s]]$f_hat_clr, col = "red")
  lines(dta_dens$y, sim_result_dat[[n]][[s]]$coverage_pw$CI_up_c, lty = 3, col = 1)
  lines(dta_dens$y, sim_result_dat[[n]][[s]]$coverage_pw$CI_low_c, lty = 3, col = 1)
  # lines(dta_dens$y, sim_result_dat[[n]][[s]]$coverage_pw$CI_up_p, lty = 3, col = 1)
  # lines(dta_dens$y, sim_result_dat[[n]][[s]]$coverage_pw$CI_low_p, lty = 3, col = 1)

  legend("bottom", bty = "n",
         legend = c(paste0("Vc: ", round(sim_result_dat[[n]][[s]]$statistic_vc, 2),
                           " (", ifelse(sim_result_dat[[n]][[s]]$coverage_Vc,
                                        "in CR", "not in CR"), ")"),
                    paste0("Vp: ", round(sim_result_dat[[n]][[s]]$statistic_vp, 2),
                           " (", ifelse(sim_result_dat[[n]][[s]]$coverage_Vp,
                                         "in CR", "not in CR"), ")"),
                    "", "truth", "estimate", "pointwise CI (Vc)"),
         lty = c(0, 0, 0, 1, 1, 3), col = c(rep("black", 3), "green", "red", "black"))

  },
  s = slider(1, S, label = "Iteration index (sorted by Vc)"),
  N = picker(5000, 10000, 50000, 100000, label = "Sample Size"),
  n_bins = picker(100, 1000, 2000, label = "Bins")
  # , type = picker("Sort by Vc", "Sort by Vp", "Unsorted") # works, but is super slow
)

# # Threshold for test decision
# K_T <- length(sim_result[[1]][[1]]$model$coefficients)
# qchisq(0.95, df = K_T - 1) # - number of covariate combinations (intercepts), here: 1

# Summary:
# The empirical coverage rates only get close to the theoretical ones with very
# large obervation numbers (>= 50000). This is also the case for the average-across-
# the-function confidence bands proposed by Wood (2017), as the plot
# "Simulation_CI_average_across_function.pdf" created above shows, and not only
# for our newly developed approach.

# Further Ideas:
# - Use clean integrate-to-zero constraint (so far, we simply use ti(..., mc = TRUE),
#   enforcing a sum-to-zero constraint, which only approximates the integrate-to-
#   zero-constraint)
# - Sample from Poisson distribution
# - Try different true density
# - Check average-across-the-function confidence bands for different gam-models
#   (not our density-setting), i.e., with global intercept and without integrate-
#   to-zero constraint, different families; Is the number of observations needed
#   still very large?
# - Create manipulate plot as above, but sorted by average-across-the-function to
#   see, which iterations are bad; Are those close to the penalty null space?
#   (First: How does the penalty null space look like in L^2_0? -> Sonja and Fabian
#   wrote a discussion paper, where they showed that Wood's confidence bands
#   perform bad, when close to penalty null space)
