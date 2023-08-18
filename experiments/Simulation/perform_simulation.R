################################################################################
################################## Functions ###################################
################################################################################

# bs_L20 creates a B-spline basis transformed to fulfill the integrate-to-zero
# constraint. Arguments are the same as in splines:splineDesign:
# - knots: a numeric vector of (inner and outer) knot positions (which will be
#   sorted increasingly if needed)
# - x: a numeric vector of values at which to evaluate the basis functions or
#   derivatives. The values in x must be between the “inner” knots knots[ord] and
#   knots[ length(knots) - (ord-1)].
# - ord: a positive integer giving the order of the spline function. This is the
#   number of coefficients in each piecewise polynomial segment, thus a cubic
#   spline has order 4. Defaults to 4.
# Value: A matrix with length(x) rows and length(knots) - ord - 1 columns. The i'th
# row of the matrix contains the evaluations of the constrained spline functions
# (defined by the knot vector and the order) at the i'th value of x.

bs_L20 <- function(knots, x, ord = 4) {
  C <- sapply(1:(length(knots) - ord), function(j) integrate(function(x)
    splines::splineDesign(knots = knots, x, ord = ord, derivs = 0)[, j],
    lower = knots[ord], upper = knots[length(knots) - (ord-1)])$value)
  Z <- MASS::Null(C)
  X_L20 <- splines::splineDesign(knots = knots, x, ord = ord, derivs = 0) %*% Z
}

# true density to interpolate (we use a beta distribution)
f_0 <- function(t, a, b) {
  dbeta(t, a, b)
}

# clr transformation of true density
f_0_clr <- function(t, a = 3, b = 3) {
  int_log_f_0 <- integrate(function(x, a, b) log(f_0(x, a, b)),
                           lower = 0, upper = 1, a = a, b = b)$value
  log(f_0(t, a, b)) - int_log_f_0
}

# interpolate_clr constructs the spline basis expansion using L^2_0 basis functions
# for a given coefficient vector theta
interpolate_clr <- function(theta, x, knots, ord = 4) {
  bs_L20(knots, x, ord = ord) %*% theta
}

# interpolate applies the inverse clr transformation of interpolate_clr
interpolate <- function(theta, x, knots, ord = 4) {
  int_exp_f_clr <- integrate(function(x, theta, knots, ord)
    exp(interpolate_clr(theta, x, knots, ord)),
    lower = 0, upper = 1, theta = theta, knots = knots, ord = ord)$value
  exp(interpolate_clr(theta, x, knots, ord)) / int_exp_f_clr
}

# get_knots computes regular knots for usage with a (constrained) B-spline basis.
# Arguments:
# - K_T: Number of B-spline functions to be used (before transformation to L^2_0)
# - ord: spline order as in splines::splineDesign (number of non-zero splines at
#        each value -> degree = ord - 1)
get_knots <- function(K_T = 10, ord = 4) {
  k_in <- seq(0, 1, length.out = K_T - 2)
  k <- c(k_in[1] - ((ord - 1):1) * diff(k_in)[1], k_in,
         k_in[length(k_in)] + (1:(ord-1)) * diff(k_in)[1])
  return(k)
}

# get_interpolation_values computes random values such that between each pair of
# knots, there is an observation.
# Arguments:
# - K_T: Number of B-spline functions to be used (before transformation to L^2_0)
# - ord: spline order as in splines::splineDes
get_interpolation_values <- function(knots, ord = 4) {
  k_in <- knots[ord:(length(knots) - ord + 1)]
  n_obs <- length(knots) - ord - 1
  t <- numeric(n_obs)
  while(any(duplicated(t))) {
    t <- sapply(seq_along(k_in[-1]), function(i) runif(1, min = k_in[i], max = k_in[i + 1]))
    t <- sort(round(c(t, runif(n_obs - length(t))), 2))
  }
  return(t)
}

# get_theta computes the coefficient vector interpolating a beta distribution.
# Arguments:
# - a, b: shape parameters of the beta distribution to be interpolated
# - K_T: Number of B-spline functions to be used (before transformation to L^2_0)
# - ord: spline order as in splines::splineDesign (number of non-zero splines at
#     each value -> degree = ord - 1)
# - knots: vector of length K_T + ord giving the inner and outer knots to be used
#     for constructing the transformed B-splines. Defaults to NULL, in which case
#     equidistant knots over [0, 1] (with equidistant outer knot extension) are used
# - t: vector of length K_T - 1 containing the values, where the function should
#     be interpolated. There should be at least one value between each adjacent
#     pair of knots. Furthermore, equidistant values seem to yield design matrices
#     that don't have full rank, causing the interpolation to fail. In this case,
#     the default (NULL) is used, which means appropriate random values are used
# Value: parameter vector, values used for interpolation, and number of trys
get_theta <- function(a = 3, b = 3, K_T = 10, ord = 4, knots = NULL, t = NULL, count = 5) {
  if (is.null(knots)) {
    knots <- get_knots(K_T, ord)
  }
  stopifnot("Number of knots must be K_T + ord!" = length(knots) == K_T + ord)

  if (is.null(t)) {
    t <- get_interpolation_values(knots = knots, ord = ord)
  }
  stopifnot("Number of values for interpolation has to be K_T - 1!" = length(t) == K_T - 1)
  t <- sort(t)

  c <- 1

  while (c < count) {
    y_clr <- f_0_clr(t, a, b)
    X <- bs_L20(knots, t, ord = ord)
    if (Matrix::rankMatrix(X)[[1]] != K_T - 1) { # Check, whether X has full rank
      warning("Resulting design matrix does not have full rank.")
    }
    theta <- try(solve(X, y_clr), silent = TRUE)
    # Check, whether theta interpolates function values. Problems seem to appear,
    # if X does not have full rank (which happens, when a regular grid for t is used,
    if (class(theta) == "try-error") {
      success <- FALSE
    } else if (!all.equal(c(X %*% theta), y_clr)) {
      success <- FALSE
    } else {
      success <- TRUE
    }
    if (!success) {
      warning("Interpolation failed, now using random values for t.")
      t <- get_interpolation_values(knots = knots, ord = ord)
      c <- c + 1
    } else {
      return(list(theta = theta, t = t, count = c))
    }
  }
  warning("Interpolation not successful.")
  return(list(theta = NA, t = t, count = c))
}

# sample_f_0_interp samples from a density constructed as inverse clr transformation
# of a constrained B-Spline basis expansion given coefficients theta by partitioning
# the support into subintervals and sampling from multinomial distribution with
# probabilities corresponding to integral; Within the bins, a uniform distribution
# is used additionally to sample values
# Arguments:
# - N: sample size
# - theta, knots, ord: vector of coefficients, vector of knots, and spline order
#     determining density to sample from
# - bins: vector containing partition of [0, 1] for the discretization of the density
sample_f_0_interp <- function(N, theta, knots, ord = 4,
                              bins = seq(knots[ord], knots[(length(knots) - (ord - 1))],
                                         length.out = 101)) {
  n_bins <- length(bins) - 1
  probs <- sapply(seq_len(n_bins), function(i)
    integrate(function(x, theta, knots, ord)
      interpolate(theta, x, knots, ord),
      theta = theta, knots = k, ord = ord,
      lower = bins[i], upper = bins[i + 1])$value)
  x <- rep(NA, N)
  for (i in seq_len(N)) {
    # x_bin <- which(rmultinom(1, 1, prob = probs) == 1)
    x_bin <- sample(n_bins, 1, prob = probs)
    x[i] <- runif(1, min = bins[x_bin], max = bins[x_bin + 1])
  }
  return(x)
}

# run_simulation performs one iteration of a simulation of the coverage rate.
# Arguments:
# - i: iteration index (for usage in lapply or even parLapply, if the simulation
#     gets computationally expensive)
# - sample_type: Either "multinomial" or "density". If "multinomial", the counts
#     for the poisson regression model are sampled directly from a multinomial
#     distribution with weights exp(eta) / sum(exp(eta)), where eta is the predictor
#     specified via step_size, theta, knots, and ord; If "density", sample_f_0_interp
#     is used to sample from the density specified via theta, knots, and ord and
#     a histogram is constructed. For stepsize -> 0 results of both approaches
#     converge to each other; Thus, for small step_size "multinomial" should be
#     used, since it is faster. Furthermore, in this case the theoretical results
#     of Wood (2017), Section 6, should hold
# - theta, knots, ord: vector of coefficients, vector of knots, and spline order
#     determining predictor (if sample_type = "multinomial") or density (if
#     sample_type = "density").
# - pen_ord: order of differences to penalize in gam() (corresponds to m[2] in gam())
# - N: sample size
# - step_size: width of histogram bins
# - alpha: determinces significance level (1 - alpha) for checking whether confidence
#     region covers true theta
run_simulation <- function(i, sample_type = c("multinomial", "density"), theta,
                           knots, ord = 4, pen_ord = 2, N = 10000, step_size = 0.01,
                           alpha = 0.05) {

  sample_type <- match.arg(sample_type)

  grid_hist <- seq(from = 0, to = 1, by = step_size)
  # results of "multinomial" converge to the ones of "density" for step_size -> 0
  if (sample_type == "multinomial") {
    mids_dens <- grid_hist[1:(length(grid_hist) - 1)] + step_size / 2
    eta <- c(log(step_size) + interpolate_clr(theta, mids_dens, knots = knots, ord = 4))
    probs <- (exp(eta)) / (sum(exp(eta)))
    counts_dens <- rmultinom(n = 1, size = N, prob = probs)
  } else {
    y <- sample_f_0_interp(N, theta, knots, ord, bins = grid_hist)
    # n_bins <- 1 / step_size
    hist_dens <- hist(y, breaks = grid_hist, right = FALSE, plot = FALSE)
    counts_dens <- hist_dens$counts
    mids_dens <- hist_dens$mids
  }

  dta_dens <- as.data.frame(cbind(counts_dens, mids_dens))
  colnames(dta_dens) <- c("counts", "y")

  K_T <- length(knots) - ord
  Delta <- rep(step_size, length(mids_dens))
  # Here, we need an intercept, which represents our one "covariate combination"
  model <- gam(counts ~ 1 + ti(y, bs = "ps", m = list(c(ord - 2, pen_ord)),
                               mc = TRUE, np = FALSE, k = K_T) + offset(log(Delta)), # for the continuous case, the sum-to-zero constraint specified via mc = TRUE is an approximation of the integrate-to-zero constraint
               data = dta_dens, knots = list(y = knots), method = "REML", family = poisson())

  X <- model.matrix(model)[, 2:K_T]
  theta_hat <- model$coefficients[2:K_T]

  f_hat_clr <- c(X %*% theta_hat)
  f_hat <- exp(f_hat_clr) / mean(exp(f_hat_clr))

  # theta ~ N(theta_hat, Vc)
  # => CR = {theta : (theta - theta_hat)^T V^{-1} (theta - theta_hat) <= chi^2_{1-alpha}(K_T)}

  # without intercepts per covariate combination (here: 1)
  Vc <- model$Vc[2:K_T, 2:K_T, drop = FALSE]
  theta_diff <- matrix(theta - theta_hat, ncol = 1)
  Vc_inv <- solve(Vc)
  chi_statistic_Vc <- t(theta_diff) %*% Vc_inv %*% theta_diff
  check_coverage_Vc <- as.numeric(chi_statistic_Vc) <= qchisq(1-alpha, df = K_T - 1) # - number of covariate combinations (intercepts), here: 1

  Vp <- model$Vp[2:K_T, 2:K_T, drop = FALSE]
  Vp_inv <- solve(Vp)
  chi_statistic_Vp <- t(theta_diff) %*% Vp_inv %*% theta_diff
  check_coverage_Vp <- as.numeric(chi_statistic_Vp) <= qchisq(1-alpha, df = K_T - 1)

  se_p <- sqrt(rowSums((X %*% Vp) * X))
  se_c <- sqrt(rowSums((X %*% Vc) * X))
  CI_up_p <- f_hat_clr + qnorm(1 - alpha / 2) * se_p
  CI_low_p <- f_hat_clr - qnorm(1 - alpha / 2) * se_p
  CI_up_c <- f_hat_clr + qnorm(1 - alpha / 2) * se_c
  CI_low_c <- f_hat_clr - qnorm(1 - alpha / 2) * se_c
  f_true_clr <- interpolate_clr(theta, dta_dens$y, knots)
  CI_p_check <- (CI_low_p <= f_true_clr) & (f_true_clr <= CI_up_p)
  CI_c_check <- (CI_low_c <= f_true_clr) & (f_true_clr <= CI_up_c)
  CIs <- data.frame(f_true_clr, CI_low_c, CI_up_c, CI_c_check, CI_low_p, CI_up_p, CI_p_check)

  print(i)
  return(list(coverage_Vc = check_coverage_Vc, coverage_Vp = check_coverage_Vp,
              statistic_vc = chi_statistic_Vc, statistic_vp = chi_statistic_Vp,
              data = dta_dens, model = model, f_hat = f_hat,
              f_hat_clr = f_hat_clr, coverage_pw = CIs))
}

################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

library(mgcv)

a <- 3
b <- 3
K_T <- 10 # 8 # 10
ord <- 4 # spline order as in splines::splineDesign (number of non-zero splines
         # at each value -> degree = ord - 1)
# m <- ord - 2 # spline order as in mgcv::gam (m = 2: cubic splines -> degree = m + 1)

set.seed(1542)
k <- get_knots(K_T = K_T, ord = ord)
t <- get_interpolation_values(knots = k, ord = ord)
theta <- get_theta(a = a, b = b, K_T = K_T, ord = ord, knots = k, t = t)$theta
# Regular grids like t <- seq(0.05, 0.95, length.out = K_T - 1)) yield design
# matrices that don't have full rank
# get_theta(a = a, b = b, K_T = K_T, ord = ord, knots = k, t = t)
# Strange: Before the transformation, the matrix did have full rank:
# test <- splines::splineDesign(knots = k, t, ord = 4, derivs = 0)
# Matrix::rankMatrix(test)

grid <- seq(0, 1, length.out = 100)
f_0_clr_interp <- interpolate_clr(theta, grid, k)
saveRDS(f_0_clr_interp, "./Simulation_Objects/f_0_clr_interp.rds")
# plot(grid, f_0_clr_interp, type = "l", ylab = "clr-density", xlab = "t",
#      main = "True clr-density")
# lines(grid, f_0_clr(grid, a = a, b = b), col = 2)
# points(t, f_0_clr(t, a, b), pch = 4)
# legend("bottom", legend = c("True clr-density (interpolation)",
#                             paste0("Beta(", a, ", ", b, ")-clr-density")),
#        col = 1:2, lty = 1, bty = "n")

# # Check:
# bs_theta <- sapply(seq_along(theta), function(i) bs_L20(k, grid)[, i] * theta[i])
# f_0_clr_interp2 <- rowSums(bs_theta) #  c(bs_L20(k, grid) %*% theta)
# matplot(grid, bs_L20(k, grid), type = "l", ylim = range(bs_theta, f_0_clr_interp2))
# matlines(grid, bs_theta)
# matlines(grid, f_0_clr_interp2)
# lines(grid, f_0_clr_interp, col = 2, lty = 2)


f_0_interp <- interpolate(theta, grid, k)
saveRDS(f_0_interp, "./Simulation_Objects/f_0_interp.rds")
# plot(grid, f_0_interp, type = "l", ylab = "density", xlab = "t",
#      main = "True density")
# lines(grid, f_0(grid, a = a, b = b), col = 2)
# legend("bottom", legend = c("True density (interpolation)",
#                             paste0("Beta(", a, ", ", b, ")-density")),
#        col = 1:2, lty = 1, bty = "n")

################################################################################
############################## Perform Simulation ##############################
################################################################################

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

alpha <- 0.05 # 95% confidence

S <- 200 # simulation runs
N_vec <- c(5000, 10000, 50000, 100000)
step_size_vec <- c(0.01, 0.001, 0.0005)
type_vec <- c("Vc", "Vp")
lengths_vecs <- sapply(list(N_vec, step_size_vec, type_vec), length)
n_bins <- sapply(step_size_vec, function(s) length(seq(0, 1, by = s))) - 1
# coverage_Vc <- coverage_Vp <- matrix(numeric(length(N_vec) * length(step_size_vec)),
#                     nrow = length(N_vec), dimnames = list(paste0("N: ", N_vec),
#                                                           paste0("Bins: ", n_bins)))

coverage_rate <- array(numeric(prod(lengths_vecs)), dim = lengths_vecs,
                       dimnames = list(paste0("N: ", N_vec),
                                       paste0("Bins: ", n_bins),
                                       paste0("Type: ", type_vec)))

set.seed(1326)
count <- 1
for(i in seq_along(N_vec)) {
  for(j in seq_along(step_size_vec)) {
    print(paste("Stepsize: ", step_size_vec[j], "; N: ", N_vec[i], "; Count: ", count))
    sim_result <- lapply(1:S, run_simulation, sample_type = "multinomial",
                         theta = theta, knots = k, N = N_vec[i],
                         step_size = step_size_vec[j], alpha = alpha)
    saveRDS(sim_result, paste0("./Simulation_Objects/S", S, "_sample_size",
                               N_vec[i], "_bins", n_bins[j], ".rds"))
    coverage_rate[i, j, 1] <- sum(sapply(1:S, function(k) sim_result[[k]]$coverage_Vc)) / S
    coverage_rate[i, j, 2] <- sum(sapply(1:S, function(k) sim_result[[k]]$coverage_Vp)) / S
    count <- count + 1
  }
}

saveRDS(coverage_rate, "./Simulation_Objects/coverage_rates.rds")
options(scipen = scipen)
