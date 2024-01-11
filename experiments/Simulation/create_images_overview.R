################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size
# options(scipen = scipen)

grid <- seq(0, 1, length.out = 100)
N_vec <- c(5000, 10000, 50000, 100000)
step_size_vec <- c(0.01, 0.001, 0.0005)
n_bins <- sapply(step_size_vec, function(s) length(seq(0, 1, by = s))) - 1
N_long <- rep(N_vec, each = length(n_bins))
n_bins_long <- rep(n_bins, length(N_vec))
S <- 200 # simulation runs
alpha <- 0.05

library(plot.matrix)
plot_results <- function(result, S = 200, main = "Coverage rates",
                         breaks = seq(0, 1, 0.1), n_bins = c(100, 1000, 2000),
                         N_vec = c(5000, 10000, 50000, 100000)) {
  plot(result, xlab = "", ylab = "Sample size N", axis.col = NULL,
       axis.row = NULL, breaks = breaks,
       main = main)
  axis(1, at = 1:3, labels = n_bins)
  axis(2, at = 4:1, labels = N_vec, las = 2)
  text(rep(1:3, each = 4), rep(1:4, 3), labels = result[nrow(result):1, ])
  mtext("Bin number G", side = 1, line = 3, at = dim(result)[2] / 2 + 0.5)
}

################################################################################
############################# Penalized estimation #############################
################################################################################

scenarios <- c("a3_b3", "a2_b5", "a2_b5_K30", "line_a01")
folder <- paste0("./penalized/", scenarios)
dist <- c("Beta(3, 3)", "Beta(2, 5)", "Beta(2, 5)", "0.95+0.1*t")
library(latex2exp)
K <- list(TeX("K_T = 9"), TeX("K_T = 9"), TeX("K_T = 29"), TeX("K_T = 9")) # c(10, 10, 30, 10)
name <- paste0("pen_", scenarios, "_")
legend_pos <- c("bottom", "topright", "topright", "bottom")
main <- list(TeX("(a) Beta(3, 3), K_T = 9"), TeX("(b) Beta(2, 5), K_T = 9"),
             TeX("(c) Beta(2, 5), K_T = 29"), TeX("(d) $0.95 + 0.1 \\cdot t$, K_T = 9"))

pdf("./Overview_Images/true_densities_pen.pdf", width = 10, height = 5)
layout(matrix(1:8, nrow = 2))
for (f in seq_along(folder)) {
  f_0_clr <- readRDS(paste0(folder[f], "/Simulation_Objects/f_0_clr_interp.rds"))
  f_0 <- readRDS(paste0(folder[f], "/Simulation_Objects/f_0_interp.rds"))
  plot(grid, f_0, type = "l", ylab = "density", xlab = "t",
       main = paste0("(", letters[f], ") Interp. of ", dist[f]), ylim = c(0, max(f_0)))
  legend(legend_pos[f], legend = K[[f]], bty = "n")
  plot(grid, f_0_clr, type = "l", ylab = "clr-density", xlab = "t",
       main = paste0("(", letters[f], ") Interp. of ", dist[f]))
  abline(h = 0, col = "grey", lty = 2)
  legend("bottom", legend = K[[f]], bty = "n")
}
dev.off()

pdf("./Overview_Images/coverage_pen.pdf", width = 13, height = 3)
layout(matrix(1:4, nrow = 1, byrow = TRUE))
par(mar = c(mar = c(5, 4.7, 4, 4) + 0.1), mgp = c(3.7, 1, 0))
for (f in seq_along(folder)) {
  coverage_rates <- readRDS(paste0(folder[f], "/Simulation_Objects/coverage_rates.rds"))
  plot_results(coverage_rates[, , 1], S = 200, main = main[[f]])
}
dev.off()

for (f in seq_along(folder)) {
  sim_result <- lapply(1:length(N_long), function(i)
    readRDS(paste0(folder[f], "/Simulation_Objects/S", S, "_N", N_long[i], "_G", n_bins_long[i], ".rds")))

  # relMSE <- readRDS(paste0(folder[f], "/Simulation_Objects/relMSE.rds"))
  relMSE <- lapply(seq_along(N_long),
                   function(i) sapply(1:S,
                                      function(j) sim_result[[i]][[j]]$relMSE))
  pdf(paste0(folder[f], "/Simulation_Images/", name[f], "relMSE.pdf"), height = 6.5, width = 5.5)
  par(mfrow = c(4, 3), mar = c(1, 4, 4, 1) + 0.1)
  sapply(1:length(N_long), function(i) {
    main <- paste0("(", letters[f], ") ", paste0("N: ", N_long[i], ", G: ", n_bins_long[i]))
    boxplot(relMSE[[i]], main = main, ylim = range(relMSE), las = 1)
    # abline(h = mean(relMSE[[i]]), col = "red")
  })
  dev.off()

  coverage_pw <- lapply(seq_along(N_long),
                        function(i) sapply(1:S,
                                           function(j) sum(sim_result[[i]][[j]]$coverage_pw$CI_c_check) / n_bins_long[i]))

  pdf(paste0(folder[f], "/Simulation_Images/", name[f], "average_across_function.pdf"), height = 6.5, width = 5.5)
  par(mfrow = c(4, 3), mar = c(1, 4, 4, 1) + 0.1)
  sapply(1:length(N_long), function(i) {
    main <- paste0("(", letters[f], ") ", paste0("N: ", N_long[i], ", G: ", n_bins_long[i]))
    boxplot(coverage_pw[[i]], main = main, las = 1, ylim = c(0, 1))
    abline(h = 1 - alpha, col = "green")
    abline(h = mean(coverage_pw[[i]]), col = "red")
  })
  dev.off()
}

################################################################################
############################ Unpenalized estimation ############################
################################################################################

scenarios <- c("a3_b3", "a2_b5", "line_a01")
folder <- paste0("./unpenalized/", scenarios)
# f_true <- list(dbeta(grid, 3, 3), dbeta(grid, 2, 5), dbeta(grid, 2, 5), 0.95 + 0.1 * grid )
dist <- c("Beta(3, 3)", "Beta(2, 5)", "0.95+0.1*t")
K <- paste0("K = ", c(10, 10, 10))
name <- paste0(scenarios, "_")
my_letters <- c("a", "b", "d")

pdf("./Overview_Images/true_densities.pdf", width = 7.3, height = 5)
layout(matrix(1:6, nrow = 2))
for (f in seq_along(folder)) {
  f_0_clr <- readRDS(paste0(folder[f], "/Simulation_Objects/f_0_clr_interp.rds"))
  f_0 <- readRDS(paste0(folder[f], "/Simulation_Objects/f_0_interp.rds"))
  plot(grid, f_0, type = "l", ylab = "density", xlab = "t",
       main = paste0("(", letters[f], ") Interp. of ", dist[f]), ylim = c(0, max(f_0)))
  legend("bottom", legend = K[f], bty = "n")
  plot(grid, f_0_clr, type = "l", ylab = "clr-density", xlab = "t",
       main = paste0("(", my_letters[f], ") Interp. of ", dist[f]))
  abline(h = 0, col = "grey", lty = 2)
  legend("bottom", legend = K[f], bty = "n")
}
dev.off()


pdf("./Overview_Images/coverage.pdf", width = 9.8, height = 3)
layout(matrix(1:3, nrow = 1, byrow = TRUE))
par(mar = c(mar = c(5, 4.7, 4, 4) + 0.1), mgp = c(3.7, 1, 0))
for (f in seq_along(folder)) {
  coverage_rates <- readRDS(paste0(folder[f], "/Simulation_Objects/coverage_rates.rds"))
  main <- paste0("(", my_letters[f], ") ", dist[f], ", ", K[f])
  plot_results(coverage_rates[, , 2], S = 200, main = main)
}
dev.off()

for (f in seq_along(folder)) {
  sim_result <- lapply(1:length(N_long), function(i)
    readRDS(paste0(folder[f], "/Simulation_Objects/S", S, "_N", N_long[i], "_G", n_bins_long[i], ".rds")))

  # relMSE <- readRDS(paste0(folder[f], "/Simulation_Objects/relMSE.rds"))
  relMSE <- lapply(seq_along(N_long),
                   function(i) sapply(1:S,
                                      function(j) sim_result[[i]][[j]]$relMSE))
  cutoff <- max(sapply(relMSE, function(r) boxplot(r)$stats[5, 1]))
  if (cutoff > 6) {
    ylim <- c(0, cutoff)
    outlier <- TRUE
  } else {
    ylim <- range(relMSE)
    outlier <- FALSE
  }

  pdf(paste0(folder[f], "/Simulation_Images/", name[f], "relMSE.pdf"), height = 6.5, width = 5.5)
  par(mfrow = c(4, 3), mar = c(1, 4, 4, 1) + 0.1)
  sapply(1:length(N_long), function(i) {
    # cutoff <- 7
    main <- paste0("(", my_letters[f], ") ", paste0("N: ", N_long[i], ", G: ", n_bins_long[i]))
    boxplot(relMSE[[i]], main = main, las = 1, ylim = ylim)
    # abline(h = mean(relMSE[[i]]), col = "red")
    if (outlier & max(relMSE[[i]]) > cutoff + 0.02 * cutoff) {
      m <- max(relMSE[[i]])
      o <- length(which(relMSE[[i]] > cutoff + 0.02 * cutoff))
      # points(1, cutoff + 0.02 * cutoff, col = "red")
      text(1, 0.85 * cutoff, paste0(o, " more values,\nmax: ", round(m)), col = 2)
    }
  })
  dev.off()

  coverage_pw <- lapply(seq_along(N_long),
                        function(i) sapply(1:S,
                                           function(j) sum(sim_result[[i]][[j]]$coverage_pw$CI_p_check) / n_bins_long[i]))

  pdf(paste0(folder[f], "/Simulation_Images/", name[f], "average_across_function.pdf"), height = 6.5, width = 5.5)
  par(mfrow = c(4, 3), mar = c(1, 4, 4, 1) + 0.1)
  sapply(1:length(N_long), function(i) {
    main <- paste0("(", my_letters[f], ") ", paste0("N: ", N_long[i], ", G: ", n_bins_long[i]))
    boxplot(coverage_pw[[i]], main = main, las = 1, ylim = c(0, 1))
    abline(h = 1 - alpha, col = "green")
    abline(h = mean(coverage_pw[[i]]), col = "red")
  })
  dev.off()
}

################################################################################
################## Compare relMSEs penalized vs. unpenalized ###################
################################################################################

scenarios <- c("a3_b3", "a2_b5", "a2_b5_K30", "line_a01",
               "a3_b3", "a2_b5", "line_a01")
folder <- c(paste0("./penalized/", scenarios[1:4]), paste0("./unpenalized/", scenarios[5:7]))
name <- paste0(scenarios[c(1, 2, 4)])

relMSE <- list()

for (f in seq_along(folder)) {
  sim_result <- lapply(1:length(N_long), function(i)
    readRDS(paste0(folder[f], "/Simulation_Objects/S", S, "_N", N_long[i], "_G", n_bins_long[i], ".rds")))

  relMSE[[f]] <- lapply(seq_along(N_long),
                   function(i) sapply(1:S,
                                      function(j) sim_result[[i]][[j]]$relMSE))
}

pairs <- list(c(1, 5), c(2, 6), c(4, 7))
for (p in seq_along(pairs)) {
  pen <- pairs[[p]][1]
  unpen <- pairs[[p]][2]
  pdf(paste0("./Overview_Images/relMSE_comparison_", name[p], ".pdf"), height = 6.5, width = 5.5)
  par(mfrow = c(4, 3), mar = c(1.9, 4, 4, 1) + 0.1)
  sapply(1:length(N_long), function(i) {
    cutoff <- max(sapply(list(relMSE[[pen]][[i]], relMSE[[unpen]][[i]]),
                         function(r) boxplot(r, plot = FALSE)$stats[5, 1]))
    ylim <- c(0, cutoff)
    main <- paste0("(", my_letters[p], ") ", paste0("N: ", N_long[i], ", G: ", n_bins_long[i]))
    boxplot(list(relMSE[[pen]][[i]], relMSE[[unpen]][[i]]), xaxt = "n", las = 1,
            main = main, ylim = ylim)
    axis(1, at = c(1, 2), labels = c("pen.", "unpen."))
  })
  dev.off()
}

################################################################################
########################## Unpenalized with smaller N ##########################
################################################################################

scenarios <- c("a3_b3", "a2_b5", "line_a01")
folder <- paste0("./unpenalized/", scenarios)
# f_true <- list(dbeta(grid, 3, 3), dbeta(grid, 2, 5), dbeta(grid, 2, 5), 0.95 + 0.1 * grid )
dist <- c("Beta(3, 3)", "Beta(2, 5)", "0.95+0.1*t")
K <- paste0("K = ", c(10, 10, 10))
name <- paste0(scenarios, "_small_N_")
my_letters <- c("a", "b", "d")
N_vec <- c(100, 500, 1000, 2500)
N_long <- rep(N_vec, each = length(n_bins))

pdf("./Overview_Images/coverage_small_N.pdf", width = 9.8, height = 3)
layout(matrix(1:3, nrow = 1, byrow = TRUE))
par(mar = c(mar = c(5, 4.7, 4, 4) + 0.1), mgp = c(3.7, 1, 0))
for (f in seq_along(folder)) {
  coverage_rates <- readRDS(paste0(folder[f], "/Simulation_Objects/coverage_rates_small_N.rds"))
  main <- paste0("(", my_letters[f], ") ", dist[f], ", ", K[f])
  plot_results(coverage_rates[, , 2], S = 200, main = main, N_vec = N_vec,
               n_bins = n_bins)
}
dev.off()

for (f in seq_along(folder)) {
  sim_result <- lapply(1:length(N_long), function(i)
    readRDS(paste0(folder[f], "/Simulation_Objects/S", S, "_N", N_long[i], "_G", n_bins_long[i], ".rds")))

  # relMSE <- readRDS(paste0(folder[f], "/Simulation_Objects/relMSE.rds"))
  relMSE <- lapply(seq_along(N_long),
                   function(i) sapply(1:S,
                                      function(j) sim_result[[i]][[j]]$relMSE))
  cutoff <- max(sapply(relMSE, function(r) boxplot(r)$stats[5, 1]))
  if (cutoff > 6) {
    ylim <- c(0, cutoff)
    outlier <- TRUE
  } else {
    ylim <- range(relMSE)
    outlier <- FALSE
  }

  pdf(paste0(folder[f], "/Simulation_Images/", name[f], "relMSE.pdf"), height = 6.5, width = 5.5)
  par(mfrow = c(4, 3), mar = c(1, 4, 4, 1) + 0.1)
  sapply(1:length(N_long), function(i) {
    # cutoff <- 7
    main <- paste0("(", my_letters[f], ") ", paste0("N: ", N_long[i], ", G: ", n_bins_long[i]))
    boxplot(relMSE[[i]], main = main, las = 1, ylim = ylim)
    # abline(h = mean(relMSE[[i]]), col = "red")
    if (outlier & max(relMSE[[i]]) > cutoff + 0.02 * cutoff) {
      m <- max(relMSE[[i]])
      o <- length(which(relMSE[[i]] > cutoff + 0.02 * cutoff))
      # points(1, cutoff + 0.02 * cutoff, col = "red")
      text(1, 0.85 * cutoff, paste0(o, " more values,\nmax: ", round(m)), col = 2)
    }
  })
  dev.off()

  coverage_pw <- lapply(seq_along(N_long),
                        function(i) sapply(1:S,
                                           function(j) sum(sim_result[[i]][[j]]$coverage_pw$CI_p_check) / n_bins_long[i]))

  pdf(paste0(folder[f], "/Simulation_Images/", name[f], "average_across_function.pdf"), height = 6.5, width = 5.5)
  par(mfrow = c(4, 3), mar = c(1, 4, 4, 1) + 0.1)
  sapply(1:length(N_long), function(i) {
    main <- paste0("(", my_letters[f], ") ", paste0("N: ", N_long[i], ", G: ", n_bins_long[i]))
    boxplot(coverage_pw[[i]], main = main, las = 1, ylim = c(0, 1))
    abline(h = 1 - alpha, col = "green")
    abline(h = mean(coverage_pw[[i]]), col = "red")
  })
  dev.off()
}
