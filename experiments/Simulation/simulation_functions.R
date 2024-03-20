################################################################################
################################## Functions ###################################
################################################################################

source("./densities.R")

### Functions to construct true underlying density in the span of a transformed
### B-Spline basis, i.e., of the form B^T theta. To obtain a nicely shaped
### true density, we interpolate a density from a beta-distribution

# get_knots computes regular knots for usage with a (constrained) B-spline basis.
# Arguments:
# - n_splines: Number of B-spline functions to be used (before transformation to L^2_0)
# - ord: spline order as in splines::splineDesign (number of non-zero splines at
#        each value -> degree = ord - 1)
get_knots <- function(n_splines = 10, ord = 4) {
  interior_knots <- seq(0, 1, length.out = n_splines - 2)
  knots <- c(interior_knots[1] - ((ord - 1):1) * diff(interior_knots)[1],
             interior_knots,
             interior_knots[length(interior_knots)] + (1:(ord-1)) * diff(interior_knots)[1])
  return(knots)
}

# get_interpolation_values computes random values such that between each pair of
# knots, there is an observation.
# Arguments:
# - n_splines: Number of B-spline functions to be used (before transformation to L^2_0)
# - ord: spline order as in splines::splineDes
get_interpolation_values <- function(knots, ord = 4) {
  interior_knots <- knots[ord:(length(knots) - ord + 1)]
  n_obs <- length(knots) - ord - 1
  quantiles <- numeric(n_obs)
  while(any(duplicated(quantiles))) {
    quantiles <- sapply(seq_along(interior_knots[-1]),
                function(i) runif(1, min = interior_knots[i],
                                  max = interior_knots[i + 1]))
    quantiles <- sort(round(c(quantiles, runif(n_obs - length(quantiles))), 2))
  }
  return(quantiles)
}


get_scenario_name <- function(density_params) {
  if (density_params$density_name == "beta") {
    scenario_name <- paste0(density_params$density_name, " ", "a=", density_params$a, " ",
                            "b=", density_params$b, " ", "n_knots=", density_params$n_knots
                            )
  }

  if (density_params$density_name == "truncated_normal") {
    scenario_name <- paste0(density_params$density_name, "_", "mean", density_params$mean, "_",
                            "sd", density_params$sd, "_", "knots", density_params$n_knots
    )
  }

  scenario_name

}

# get_theta computes the coefficient vector interpolating a beta distribution.
# Arguments:
# - a, b: shape parameters of the beta distribution to be interpolated
# - n_splines: Number of B-spline functions to be used (before transformation to L^2_0)
# - ord: spline order as in splines::splineDesign (number of non-zero splines at
#     each value -> degree = ord - 1)
# - knots: vector of length n_splines + ord giving the inner and outer knots to be used
#     for constructing the transformed B-splines. Defaults to NULL, in which case
#     equidistant knots over [0, 1] (with equidistant outer knot extension) are used
# - t: vector of length n_splines - 1 containing the values, where the function should
#     be interpolated. There should be at least one value between each adjacent
#     pair of knots. Furthermore, equidistant values seem to yield design matrices
#     that don't have full rank, causing the interpolation to fail. In this case,
#     the default (NULL) is used, which means appropriate random values are used
# Value: parameter vector, values used for interpolation, and number of trys
get_theta <- function(clr_density, n_splines = 10, ord = 4, knots = NULL, quantiles = NULL, count = 5) {
  if (is.null(knots)) {
    knots <- get_knots(n_splines = n_splines, ord = ord)
  }
  stopifnot("Number of knots must be n_splines + ord!" = length(knots) == n_splines + ord)

  if (is.null(quantiles)) {
    quantiles <- get_interpolation_values(knots = knots, ord = ord)
  }
  stopifnot("Number of values for interpolation has to be n_splines - 1!" = length(quantiles) == n_splines - 1)
  quantiles <- sort(quantiles)

  c <- 1

  while (c < count) {
    y_clr <- clr_density(quantiles)
    X <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = ord)
    if (Matrix::rankMatrix(X)[[1]] != n_splines - 1) { # Check, whether X has full rank
      warning("Resulting design matrix does not have full rank.")
    }
    theta <- try(solve(X, y_clr), silent = TRUE)
    # Check, whether theta interpolates function values. Problems seem to appear,
    # if X does not have full rank (which happens, when a regular grid for t is used)
    if (class(theta) == "try-error") {
      success <- FALSE
    # see https://stackoverflow.com/questions/18025797/invalid-argument-type-error-with-all-equal-r
    } else if (!isTRUE(all.equal(c(X %*% theta), y_clr))) {
      success <- FALSE
    } else {
      success <- TRUE
    }
    if (!success) {
      warning("Interpolation failed, now using random values for t.")
      quantiles <- get_interpolation_values(knots = knots, ord = ord)
      c <- c + 1
    } else {
      return(list(theta = theta, quantiles = quantiles, count = c))
    }
  }
  warning("Interpolation not successful.")
  return(list(theta = NA, quantiles = quantiles, count = c))
}

### Functions to perform a simulation given a true density via a B-spline basis
### and corresponding coefficients

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


# load self-written smoother "ms" for mixed reference measure
source("../mixed_density_smooth.R")
# source("../Eva_SOEP/help_functions.R")

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
# - pen_ord: order of differences to penalize in gam() (corresponds to m[2] in gam(..., ti()))
# - N: sample size
# - step_size: width of histogram bins
# - alpha: determinces significance level (1 - alpha) for checking whether confidence
#     region covers true theta
# - sp: smoothing parameter used for ti() in gam()
# - norm_true: norm of true density of simulation scenario for computing the
#     relMSE; defaults to NULL, in which case it is computed by the function
#     (since we perform simulations based on the same true density with different
#     observation numbers and bin width, the same norm would be computed several
#     times, thus, it saves time to compute it once and then pass it to the function)
# - bs: type of basis functions used for the model; Default "md" corresponds to
#     our mixed density smoother, "ad" to adaptive smooths (based on P-splines
#     where the smoothing parameter itself is evaluated smoothly via P-splines);
#     "ad" is just a "quick-and-dirty"-solution, using a sum-to-zero constraint
#     (which is an approximation of the actual integrate-to-zero constraint) to
#     check, whether this may be a solution for the penalized estimation. First
#     results indicate that it seems not to improve the coverage rates...
# - ad_m: Argument m for bs = "ad", specifying the dimension of the P-spline basis
#     used for the smoothing parameter. If bs = "md" this argument is ignored
run_simulation <- function(i, sample_type = c("multinomial", "density"), theta,
                           knots, ord = 4, pen_ord = 2, N = 10000, step_size = 0.01,
                           alpha = 0.05, sp = NULL, norm_true = NULL,
                           bs = c("md", "ad"), ad_m = 5) {

  sample_type <- match.arg(sample_type)
  bs <- match.arg(bs)

  n_splines <- length(knots) - ord
  grid_hist <- seq(from = 0, to = 1, by = step_size)
  success <- FALSE

  unpenalized <- ifelse(is.null(sp), FALSE, ifelse(sp == 0, TRUE, FALSE))

  while(!success) {
    # results of "multinomial" converge to the ones of "density" for step_size -> 0
    if (unpenalized) {
      support_check <- FALSE
      while (!support_check) {
        if (sample_type == "multinomial") {
          mids_dens <- grid_hist[1:(length(grid_hist) - 1)] + step_size / 2
          eta <- c(log(step_size) + spline_clr_density(mids_dens, theta, knots = knots, ord = 4))
          probs <- (exp(eta)) / (sum(exp(eta)))
          counts_dens <- rmultinom(n = 1, size = N, prob = probs)
        } else {
          y <- sample_f_0_interp(N, theta, knots, ord, bins = grid_hist)
          # n_bins <- 1 / step_size
          hist_dens <- hist(y, breaks = grid_hist, right = FALSE, plot = FALSE)
          counts_dens <- hist_dens$counts
          mids_dens <- hist_dens$mids
        }
        check <- sapply(seq_len(n_splines),
                        function(k) sum(counts_dens[which(mids_dens >= knots[k] &
                                                            mids_dens <= knots[k + ord])]))
        support_check <- all(check > 0)
      }
    } else {
      if (sample_type == "multinomial") {
        mids_dens <- grid_hist[1:(length(grid_hist) - 1)] + step_size / 2
        eta <- c(log(step_size) + interpolate_clr(mids_dens, theta, knots = knots, ord = 4))
        probs <- (exp(eta)) / (sum(exp(eta)))
        counts_dens <- rmultinom(n = 1, size = N, prob = probs)
      } else {
        y <- sample_f_0_interp(N, theta, knots, ord, bins = grid_hist)
        # n_bins <- 1 / step_size
        hist_dens <- hist(y, breaks = grid_hist, right = FALSE, plot = FALSE)
        counts_dens <- hist_dens$counts
        mids_dens <- hist_dens$mids
      }
    }


    dta_dens <- as.data.frame(cbind(counts_dens, mids_dens))
    colnames(dta_dens) <- c("counts", "y")

    Delta <- rep(step_size, length(mids_dens))
    xt_c <- list(values_discrete = FALSE, domain_continuous = c(0, 1))
    # Here, we need an intercept, which represents our one "covariate combination"
    if (bs == "md") {
      model <- gam(counts ~ 1 +
                     ti(y, bs = "md", m = list(c(ord - 2, pen_ord)), mc = FALSE,
                        np = FALSE, k = n_splines, xt = list(xt_c), sp = sp) +
                     offset(log(Delta)),
                   data = dta_dens, knots = list(y = knots), method = "REML", family = poisson())
    } else if (bs == "ad") {
      model <- gam(counts ~ 1 +
                     s(y, bs = "ad", m = ad_m, k = n_splines) +
                     offset(log(Delta)),
                   data = dta_dens, knots = list(y = knots), method = "REML", family = poisson())
    }

    # remove intercepts per covariate combination (here no covariates, i.e., one intercept)
    X <- model.matrix(model)[, 2:n_splines]
    theta_hat <- model$coefficients[2:n_splines]

    theta_diff <- theta - theta_hat
    MSE <- integrate(function(x) interpolate_clr(x, theta = theta_diff, knots = knots)^2, lower = 0, upper = 1)$value
    if (is.null(norm_true)) {
      norm_true <- integrate(function(x) interpolate_clr(x, theta = theta, knots = knots)^2, lower = 0, upper = 1)$value
    }
    relMSE <- MSE / norm_true

    # theta ~ N(theta_hat, Vc)
    # => CR = {theta : (theta - theta_hat)^T V^{-1} (theta - theta_hat) <= chi^2_{1-alpha}(n_splines)}

    theta_diff <- matrix(theta_diff, ncol = 1)
    # without intercepts per covariate combination (here: 1)
    if (is.null(sp)) {
      Vc <- model$Vc[2:n_splines, 2:n_splines, drop = FALSE]
      Vc_inv <- try(solve(Vc), silent = TRUE)
    } else {
      Vc_inv <-  NA
    }

    Vp <- model$Vp[2:n_splines, 2:n_splines, drop = FALSE]
    Vp_inv <- try(solve(Vp), silent = TRUE)
    success <- ifelse("try-error" %in% union(class(Vc_inv), class(Vp_inv)), FALSE, TRUE)
  }
  if (!any(is.na(Vc_inv))) {
    chi_statistic_Vc <- t(theta_diff) %*% Vc_inv %*% theta_diff
  } else {
    chi_statistic_Vc <- NA
  }

  check_coverage_Vc <- as.numeric(chi_statistic_Vc) <= qchisq(1-alpha, df = n_splines - 1) # - number of covariate combinations (intercepts), here: 1

  chi_statistic_Vp <- t(theta_diff) %*% Vp_inv %*% theta_diff
  check_coverage_Vp <- as.numeric(chi_statistic_Vp) <= qchisq(1-alpha, df = n_splines - 1)

  f_hat_clr <- c(X %*% theta_hat)
  f_hat <- try(interpolate(dta_dens$y, theta_hat, knots, ord), silent = TRUE) # FDboost::clr(f_hat_clr, w = step_size, inverse = TRUE) # exp(f_hat_clr) / mean(exp(f_hat_clr))
  if ("try-error" %in% class(f_hat)) {
    f_hat <- FDboost::clr(f_hat_clr, w = step_size, inverse = TRUE)
  } else {
    f_hat <- c(f_hat)
  }

  se_p <- sqrt(rowSums((X %*% Vp) * X))
  if (is.null(sp)) {
    se_c <- sqrt(rowSums((X %*% Vc) * X))
  } else {
    se_c <- NA
  }
  CI_up_p <- f_hat_clr + qnorm(1 - alpha / 2) * se_p
  CI_low_p <- f_hat_clr - qnorm(1 - alpha / 2) * se_p
  CI_up_c <- f_hat_clr + qnorm(1 - alpha / 2) * se_c
  CI_low_c <- f_hat_clr - qnorm(1 - alpha / 2) * se_c
  f_true_clr <- interpolate_clr(dta_dens$y, theta, knots)
  CI_p_check <- (CI_low_p <= f_true_clr) & (f_true_clr <= CI_up_p)
  CI_c_check <- (CI_low_c <= f_true_clr) & (f_true_clr <= CI_up_c)
  CIs <- data.frame(f_true_clr, CI_low_c, CI_up_c, CI_c_check, CI_low_p, CI_up_p, CI_p_check)

  print(i)
  return(list(data = dta_dens, model = model, f_hat = f_hat, f_hat_clr = f_hat_clr,
              MSE = MSE, relMSE = relMSE, coverage_Vc = check_coverage_Vc,
              coverage_Vp = check_coverage_Vp, statistic_vc = chi_statistic_Vc,
              statistic_vp = chi_statistic_Vp, coverage_pw = CIs))
}


prepare_simulation <- function(scenario_number,
                               density,
                               seeds,
                               knots_for_scenarios,
                               spline_order,
                               number_of_b_splines,
                               grid,
                               save_path) {
  set.seed(seeds[scenario_number])

  interpolation_values <- get_interpolation_values(knots = knots_for_scenarios[[scenario_number]], ord = spline_order)
  theta <- get_theta(a = a[scenario_number], b = b[scenario_number], n_splines = number_of_b_splines[scenario_number], # and here
                     ord = spline_order, knots = knots_for_scenarios[[scenario_number]], t = interpolation_values)$theta
  saveRDS(theta, paste0(save_path, "/theta.rds"))# clr transformation of true density
  f_0_clr_interp <- interpolate_clr(x = grid, theta = theta, knots = knots_for_scenarios[[scenario_number]])
  saveRDS(f_0_clr_interp, paste0(save_path, "/f_0_clr_interp.rds"))

  f_0_interp <- interpolate(x = grid, theta = theta, knots = knots_for_scenarios[[scenario_number]])
  saveRDS(f_0_interp, paste0(save_path, "/f_0_interp.rds"))

  #plot_interpolated_density(grid = grid, interpolated_density_results = f_0_interp,
  #                         true_density = f_0, scenario_number = scenario_number,
  #
}


plot_interpolated_density <- function(grid, # and here
                                      interpolated_density_results,
                                      true_density,
                                      scenario_number,
                                      a,
                                      b,
                                      interpolation_values,
                                      is_clr = F+ALSE) {
  clr_str = ""
  if (is_clr) {
    clr_str = "clr-"
  }
  par(mfrow = c(1, 2))
  plot(grid, interpolated_density_results, type = "l", ylab = paste0(clr_str, "density"), xlab = "t",
       main = paste0("True ", clr_str, "density"),
       ylim = range(interpolated_density_results, true_density(grid, a = a, b = b), finite = TRUE))
  lines(grid, true_density(grid, a = a, b = b), col = "grey")
  points(interpolation_values, true_density(interpolation_values, a, b), pch = 4)
  legend("bottom", legend = c(paste0("True ",  clr_str ,"density (interpolation)"),
                              paste0("Beta(", a, ", ", b, paste0(")-",clr_str, "density"))), # adjust for different density
         col = 1:2, lty = 1, bty = "n")
}






