################################################################################
################################## Functions ###################################
################################################################################

source("./densities.R")
library(data.table)
library(mgcv)
library("r2r")

### Functions to construct true underlying density in the span of a transformed
### B-Spline basis, i.e., of the form B^T theta. To obtain a nicely shaped
### true density, we interpolate a density from a beta-distribution

# get_knots computes regular knots for usage with a (constrained) B-spline basis.
# Arguments:
# - n_splines: Number of B-spline functions to be used (before transformation to L^2_0)
# - ord: spline order as in splines::splineDesign (number of non-zero splines at
#        each value -> degree = ord - 1)
get_knots <- function(n_splines = 10, ord = 4, range_ = c(0,1)) {
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
get_interpolation_values <- function(knots, ord = 4, n_obs = 0) {
  interior_knots <- knots[ord:(length(knots) - ord + 1)]
  if (n_obs == 0) {
    n_obs <- length(knots) - ord - 1
  }
  quantiles <- numeric(n_obs)
  # we have to avoud that the first quantile is zero because of log trafo
  while(any(duplicated(quantiles)) & any(quantiles == 0)) {
    quantiles <- sapply(seq_along(interior_knots[-1]),
                function(i) runif(1, min = interior_knots[i],
                                  max = interior_knots[i + 1]))
    quantiles <- sort(round(c(quantiles, runif(n_obs - length(quantiles))), 4))
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
# get_theta <- function(clr_density, n_splines = 10, ord = 4, knots = NULL, quantiles = NULL, count = 5) {
#   if (is.null(knots)) {
#     knots <- get_knots(n_splines = n_splines, ord = ord)
#   }
#   stopifnot("Number of knots must be n_splines + ord!" = length(knots) == n_splines + ord)
#
#   if (is.null(quantiles)) {
#     quantiles <- get_interpolation_values(knots = knots, ord = ord)
#   }
#   stopifnot("Number of values for interpolation has to be n_splines - 1!" = length(quantiles) == n_splines - 1)
#   quantiles <- sort(quantiles)
#
#   c <- 1
#
#   while (c < count) {
#     y_clr <- clr_density(quantiles)
#     X <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = ord)
#     if (Matrix::rankMatrix(X)[[1]] != n_splines - 1) { # Check, whether X has full rank
#       warning("Resulting design matrix does not have full rank.")
#     }
#     theta <- try(solve(X, y_clr), silent = TRUE)
#     # Check, whether theta interpolates function values. Problems seem to appear,
#     # if X does not have full rank (which happens, when a regular grid for t is used)
#     if (class(theta) == "try-error") {
#       success <- FALSE
#     # see https://stackoverflow.com/questions/18025797/invalid-argument-type-error-with-all-equal-r
#     } else if (!isTRUE(all.equal(c(X %*% theta), y_clr))) {
#       success <- FALSE
#     } else {
#       success <- TRUE
#     }
#     if (!success) {
#       warning("Interpolation failed, now using random values for t.")
#       quantiles <- get_interpolation_values(knots = knots, ord = ord)
#       c <- c + 1
#     } else {
#       return(list(theta = theta, quantiles = quantiles, count = c))
#     }
#   }
#   warning("Interpolation not successful.")
#   return(list(theta = NA, quantiles = quantiles, count = c))
# }

get_approx_results <-  function(clr_density, n_splines = 10, ord = 4, knots = NULL, quantiles = NULL) {
  if (is.null(knots)) {
    knots <- get_knots(n_splines = n_splines, ord = ord)
  }
  stopifnot("Number of knots must be n_splines + ord!" = length(knots) == n_splines + ord)

  if (is.null(quantiles)) {
    quantiles <- get_interpolation_values(knots = knots, ord = ord)
  }
  stopifnot("Number of values for interpolation has to be n_splines - 1!" = length(quantiles) == n_splines - 1)
  quantiles <- sort(quantiles)

  y_clr <- clr_density(quantiles)
  dat <- data.frame(y_clr)
  dat$X <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = ord)
  penalty_matrix <- diag(length(quantiles))
  print("approximating")
  model <- gam(y_clr ~ X - 1, data = dat, paraPen = list(X=list(penalty_matrix)))
  density_params <- list(theta = model$coefficients,
                         knots = knots,
                         order = 4,
                         density_name = "spline")
  spline_densities <- get_densities(density_params, calculate_norm = TRUE)
  list(theta = model$coefficients,
       knots = knots,
       spline_densities = spline_densities)
}


get_approx_results_with_covariates <-  function(density_params,
                                                covariates,
                                                n_splines = 10,
                                                ord = 4,
                                                knots_density_ = NULL,
                                                knots_smooth_covariate = NULL,
                                                quantiles_density = NULL) {
  if (is.null(knots_density_)) {
    knots_density_ <- get_knots(n_splines = n_splines, ord = ord)
  }
  if (is.null(quantiles_density)) {
    quantiles_density <- get_interpolation_values(knots = knots_density_, n_obs = 100, ord = ord)
  }

  true_densities <- get_densities_with_covariates(density_params,
                                                covariates,
                                                calculate_norm = FALSE)
  # TODO duplicated code
  if (("rep" %in% colnames(covariates))) {
    covariates <- data.table(covariates)[, rep := NULL]
  } else {
    covariates <- data.table(covariates)
  }
  unique_covariates <- covariates[,.N, by = names(covariates)]
  n_unique_cov_combis <- nrow(unique_covariates)
  dta_dens <- data.frame(matrix(ncol = 2 + ncol(covariates),
                                nrow = 0))
  colnames(dta_dens) <- c("y_clr", "quantiles", colnames(covariates))

  # Now I will blockwise construct the design matrix
  for (i in 1:nrow(covariates)) {
    y_clr <- true_densities[[i]]$clr_density_function(quantiles_density)
    covariates_for_observation <- covariates[replicate(length(y_clr), i), ]
    dta_dens_obs <- as.data.frame(cbind(y_clr, quantiles_density, covariates_for_observation))
    colnames(dta_dens_obs) <- c("y_clr", "quantiles",
                                colnames(unique_covariates)[-length(colnames(unique_covariates))])
    dta_dens <- rbind(dta_dens, dta_dens_obs)
  }

  #design_info <- get_design_matrix_and_basis_ranges_with_covariates(dta_dens, knots_density_)
  xt_c = c(0, 1)
  # estimating --> TODO
  print("approximating")
  # maybe we can delete the penalization in x direction here as well
  model <- gam(y_clr ~ -1
               + ti(quantiles, bs = "d", m = list(c(2, 2)), mc = FALSE,
                      np = FALSE, k = n_splines, xt = list(list(xt_c)))
               + ti(quantiles, bs = "d", m = list(c(2, 2)), mc = FALSE,
                    np = FALSE, k = n_splines, by = binary_variable, xt = list(list(xt_c)))
               + ti(quantiles, bs = "d", m = list(c(2, 2)), mc = FALSE,
                                    np = FALSE, k = n_splines, by = linear_variable,
                    xt = list(list(xt_c)))
               + ti(quantiles, smooth_variable, bs = c("d","ps"), m = list(c(2, 2), c(2, 2)),
                    k = c(n_splines, 8), mc = c(FALSE, TRUE), np = FALSE, xt = list(list(xt_c))),
               data = dta_dens, knots = list(quantiles = knots_density_), method = "REML")
  # TODO plotting

  list(theta = model$coefficients,
       knots = knots_density_,
       base_range = c(1, 9),
       binary_range = c(10, 18),
       linear_range = c(19, 27),
       smooth_range = c(28, length(model$coefficients)),
       knots_smooth_covariate = model$smooth[[4]]$margin[[2]]$knots)

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
#source("../mixed_density_smooth.R")
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
#     converge to each otbher; Thus, for small step_size "multinomial" should be
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
run_simulation <- function(i, sample_type = c("bin", "value"), approx_results,
                           n_knots, ord = 4, pen_ord = 2, N = 10000, step_size = 0.01,
                           alpha = 0.05, sp = NULL, norm_true = NULL,
                           bs = c("md", "ad"), ad_m = 5) {

  # brauche ich n_knots überhaupt irgendwo?
  sample_type <- match.arg(sample_type)
  bs <- match.arg(bs)
  knots <- approx_results$knots
  theta <- approx_results$theta
  spline_densities <- approx_results$spline_densities


  n_splines <- length(knots) - ord
  success <- FALSE

  unpenalized <- ifelse(is.null(sp), FALSE, ifelse(sp == 0, TRUE, FALSE))

  while(!success) {
    # results of "multinomial" converge to the ones of "density" for step_size -> 0
    density_data <- get_density_data(spline_densities, unpenalized, step_size,
                                     n_samples = N , sample_mode = "bin",
                                     n_splines = n_splines, knots = knots,
                                     order = ord)
    xt_c <- list(values_discrete = FALSE, domain_continuous = c(0, 1))
    # Here, we need an intercept, which represents our one "covariate combination"
    if (bs == "md") {
      model <- gam(counts ~ 1 +
                     ti(y, bs = "md", m = list(c(ord - 2, pen_ord)), mc = FALSE,
                        np = FALSE, k = n_splines, xt = list(xt_c), sp = sp) +
                     offset(log(density_data$Delta)),
                   data = density_data$df, knots = list(y = knots), method = "REML", family = poisson())
    } else if (bs == "ad") {
      model <- gam(counts ~ 1 +
                     s(y, bs = "ad", m = ad_m, k = n_splines) +
                     offset(log(density_data$Delta)),
                   data = density_data$df, knots = list(y = knots), method = "REML", family = poisson())
    }

    # remove intercepts per covariate combination (here no covariates, i.e., one intercept)
    X <- model.matrix(model)[, 2:n_splines]
    theta_hat <- model$coefficients[2:n_splines]
    theta_diff <- theta - theta_hat

    estimated_density_params <- list(theta = theta_hat,
                                     knots = knots,
                                     order = ord,
                                     density_name = "spline")
    estimated_spline_densities <- get_densities(estimated_density_params,
                                               calculate_norm = TRUE)

    diff_density_params <- list(theta = theta_diff,
                                knots = knots,
                                order = ord,
                                density_name = "spline")
    diff_spline_densities <- get_densities(estimated_density_params,
                                               calculate_norm = TRUE)

    MSE <- diff_spline_densities$norm_clr_density
    if (is.null(norm_true)) {
      norm_true <- spline_densities$norm_clr_density
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

  f_hat_clr <- estimated_spline_densities$clr_density_function(density_data$df$y) # möchte design matrix mitgeben können
  f_hat <- estimated_spline_densities$density_function(density_data$df$y)

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
  f_true_clr <- spline_densities$clr_density_function(density_data$df$y)
  CI_p_check <- (CI_low_p <= f_true_clr) & (f_true_clr <= CI_up_p)
  CI_c_check <- (CI_low_c <= f_true_clr) & (f_true_clr <= CI_up_c)
  CIs <- data.frame(f_true_clr, CI_low_c, CI_up_c, CI_c_check, CI_low_p, CI_up_p, CI_p_check)

  print(i)
  return(list(data = df, model = model, f_hat = f_hat, f_hat_clr = f_hat_clr,
              MSE = MSE, relMSE = relMSE, coverage_Vc = check_coverage_Vc,
              coverage_Vp = check_coverage_Vp, statistic_vc = chi_statistic_Vc,
              statistic_vp = chi_statistic_Vp, coverage_pw = CIs))
}

run_simulation_with_covariates <- function(i, sample_type = c("bin", "value"), covariates, approx_results,
                           n_knots, ord = 4, pen_ord = 2, n_obs = 10000, step_size = 0.01,
                           alpha = 0.05, sp = NULL, norm_true = NULL,
                           bs = c("md", "ad"), ad_m = 5) {

  # brauche ich n_knots überhaupt irgendwo?
  sample_type <- match.arg(sample_type)
  bs <- match.arg(bs)
  knots <- approx_results$knots
  theta <- approx_results$theta # complete theta matrix
  base_range <- approx_results$base_range
  binary_range <- approx_results$binary_range
  linear_range <- approx_results$linear_range
  smooth_range <- approx_results$smooth_range
  knots_smooth_covariate = approx_results$knots_smooth_covariate
  density_params <- list(theta = theta,
                         knots = knots,
                         knots_smooth_covariate = knots_smooth_covariate,
                         order = ord,
                         base_range = base_range,
                         binary_range = binary_range,
                         linear_range = linear_range,
                         smooth_range = smooth_range,
                         density_name = "spline")

  spline_densities <- get_densities_with_covariates(density_params = density_params,
                                             covariates = covariates
                                             )

  n_splines <- length(knots) - ord
  # how to get spline densities with more thetas (for every spline component)?
  success <- FALSE

  unpenalized <- FALSE # penalization necessary for regression

  # results of "multinomial" converge to the ones of "density" for step_size -> 0
  density_data <- get_density_data_with_covariates(spline_densities, unpenalized, step_size,
                                   covariates = covariates, sample_mode = "bin",
                                   knots = knots,
                                   order = ord)
  # Here, we need one intercept for each unique covariate combination
  model <- gam(counts ~ 1
               + ti(y, bs = "d", m = list(c(ord - 2, pen_ord)), mc = FALSE,
                    np = FALSE, k = n_splines, xt = list(c(0, 1)), sp = sp)
               + ti(y, bs = "d", m = list(c(ord - 2, pen_ord)), mc = FALSE,
                    np = FALSE, k = n_splines, xt = list(c(0, 1)), sp = sp, by = binary_variable)
               + ti(y, bs = "d", m = list(c(ord - 2, pen_ord)), mc = FALSE,
                    np = FALSE, k = n_splines, xt = list(xt_c), sp = sp, by = linear_variable)
               + ti(share, smooth_variable, bs = c("d","ps"), m = list(c(ord - 2, pen_ord), c(ord - 2, pen_ord)),
                    k = c(12, 8), mc = c(FALSE, TRUE), np = FALSE)
               + as.factor(group_id) - 1
               + offset(log(density_data$Delta)),
               data = density_data$df, knots = list(y = knots), method = "REML", family = poisson())


    # remove intercepts per covariate combination (here no covariates, i.e., one intercept)
    n_groups <- max(density_data$group_id)
    n_params <- length(model$coefficients)
    X <- model.matrix(model)[, (n_groups + 1):n_params]
    knots_smooth_covariate_estimate <- model$smooth[[4]]$margin[[2]]$knots
    theta_hat <- model$coefficients[(n_groups + 1):n_params]
    theta_diff <- theta - theta_hat

    ## TODO

    estimated_density_params <- list(theta = theta_hat,
                                     knots = knots,
                                     knots_smooth_covariate = knots_smooth_covariate_estimate,
                                     order = ord,
                                     base_range = base_range,
                                     binary_range = binary_range,
                                     linear_range = linear_range,
                                     smooth_range = smooth_range,
                                     density_name = "spline")
    estimated_spline_densities <- get_densities_with_covariates(density_params = estimated_density_params,
                                                                covariates = covariates,
                                                                calculate_norm = TRUE
    )

    diff_density_params <- list(theta = theta_diff,
                                knots = knots,
                                knots_smooth_covariate = knots_smooth_covariate_estimate,
                                order = ord,
                                base_range = base_range,
                                binary_range = binary_range,
                                linear_range = linear_range,
                                smooth_range = smooth_range,
                                density_name = "spline")
    diff_spline_densities <- get_densities_with_covariates(density_params = diff_density_params,
                                                           covariates = covariates,
                                                           calculate_norm = TRUE
    )
    # here I have continue and calculate the Overall MSE and then break down, also calc mse for each component density
    # the coverage rates for all components and the overall densitiy
    MSE <- calculate_mse(spline_densities, diff_spline_densities, norm_true)

    # theta ~ N(theta_hat, Vc)
    # => CR = {theta : (theta - theta_hat)^T V ^{-1} (theta - theta_hat) <= chi^2_{1-alpha}(n_splines)}

    theta_diff <- matrix(theta_diff, ncol = 1) # ??

    # base component
    coverage_base <- get_coverage(model, theta_diff, "base", density_params$base_range, density_params$base_range,
                                  covariates, n_groups, sp, estimated_spline_densities, n_splines, X, density_data$quantiles)
    # binary component
    coverage_binary <- get_coverage(model, theta_diff, "binary", density_params$base_range, density_params$base_range,
                                    covariates, n_groups, sp, estimated_spline_densities, n_splines, X, density_data$quantiles)
    # linear component
    coverage_linear <- get_coverage(model, theta_diff, "linear", density_params$base_range, density_params$base_range,
                                    covariates, n_groups, sp, estimated_spline_densities, n_splines, X, density_data$quantiles)
    # smooth component
    coverage_smooth <- get_coverage(model, theta_diff, "smooth", density_params$base_range, density_params$base_range,
                                    covariates, n_groups, sp, estimated_spline_densities, n_splines, X, density_data$quantiles)
    # whole density
    coverage_whole_density <- get_coverage_density(model, theta_diff, covariates, n_groups, sp, estimated_spline_densities, n_splines, X)

  print(i)
  return(list(data = df, model = model, f_hat = f_hat, f_hat_clr = f_hat_clr,
              MSE = MSE, relMSE = relMSE, coverage_Vc = check_coverage_Vc,
              coverage_Vp = check_coverage_Vp, statistic_vc = chi_statistic_Vc,
              statistic_vp = chi_statistic_Vp, coverage_pw = CIs))
}


get_coverage <- function(model, theta_diff, effect_type, param_range, base_range,
                         covariates, n_groups, sp, estimated_spline_densities, n_splines, X,
                         quantiles) {
  alpha <- 0.05
  # I don't need the range of X thetas since we only need the thetas corresponding to the y direction
  unique_covariates_for_effect <- get_unique_covariates_for_effect(covariates, effect_type) # TODO
  X_for_effect <- X[, param_range[1]:param_range[2]]
  # do I need the basis with constraint or without?
  unity_matrix <- diag(1, nrow = n_splines - 1)
  # TODO check whether this produces the basis fo reach observation individually
  basis_effect <- get_basis_for_effect(effect_type, unique_covariates_for_effect) # I want the coverage for all unique covariates
  basis_functional_intercept <- X[1, base_range[0]:base_range[1]]
  S <- get_subsetting_matrix_S(param_range, length(theta_diff))

  # for the individual effects I don't need the S_j matrix directly as I can just subset V_c or V_p respectively
  # if I want to construct the KIs for the whole density I will have to construct it explicitly
  # therefore I construct it here then I can use it also in the KIs of the individual effects (not only for thetas)
  # without intercepts per covariate combination (here: 1)
  # I'm only allowed to include theta in y-direction
  if (is.null(sp)) {
    Vc <- model$Vc[(n_groups + 1 + param_range[1]):param_range[2], (n_groups + 1 + param_range[1]):param_range[2], drop = FALSE]
    Vc_inv <- try(solve(Vc), silent = TRUE)
  } else {
    Vc_inv <-  NA
  }

  Vp <- model$Vp[(n_groups + 1 + param_range[1]):param_range[2], (n_groups + 1 + param_range[1]):param_range[2], drop = FALSE]
  Vp_inv <- try(solve(Vp), silent = TRUE)
  success <- ifelse("try-error" %in% union(class(Vc_inv), class(Vp_inv)), FALSE, TRUE)

  check_coverage_Vc <- rep(NA, length(basis_effect))
  chi_statistic_Vc <- rep(NA, length(basis_effect))
  chi_statistic_Vp <- rep(NA, length(basis_effect))
  check_coverage_Vp <- rep(NA, length(basis_effect))
  A_for_effect_bases <- hashmap()
  ki_infos <- list()

  for (i in 1:length(basis_effect)) {
    mixed_basis <- kronecker(t(basis_effect[i]), unity_matrix)
    A_for_effect_base <- mixed_basis %*% S
    theta_diff_A <- A_for_effect_base %*% theta_diff
    # I store the different As in hash table to retrieve them fast for overall density
    # KIs and coverage
    A_for_effect_bases[[unique_covariates_for_effect[i]]] <- A_for_effect_base
    Vp_mixed <- mixed_basis %*% Vp_inv %*% t(mixed_basis)
    Vp_mixed_inv <- try(solve(Vp), silent = TRUE)
    success <- ifelse("try-error" %in% union(class(Vc_inv), class(Vp_inv)), FALSE, TRUE)

    if (!any(is.na(Vc_inv))) {
      Vc_mixed <- mixed_basis %*% Vc_inv %*% t(mixed_basis)
      Vc_mixed_inv <- try(solve(Vc_mixed), silent = TRUE)
      chi_statistic_Vc[i] <- t(theta_diff_A) %*% Vc_mixed_inv %*% theta_diff_A
    } else {
      chi_statistic_Vc[i] <- NA
    }

    check_coverage_Vc[i] <- as.numeric(chi_statistic_Vc) <= qchisq(1 - alpha, df = n_splines - 1)
    chi_statistic_Vp[i] <- t(theta_diff_A) %*% Vp_mixed_inv %*% theta_diff_A
    check_coverage_Vp[i] <- as.numeric(chi_statistic_Vp) <= qchisq(1 - alpha, df = n_splines - 1)
    # calcualte KIs
    ki_infos[i] <- get_kis(estimated_spline_densities, basis_functional_intercept,
                           Vc_mixed, Vp_mixed, quantiles)

  }


  return(ki_info = ki_infos, check_coverage_Vc = check_coverage_Vc, check_coverage_Vp = check_coverage_Vp, S = S, A = A_for_effect_bases)
}


get_unique_covariates_for_effect <- function(covariates, effect_type) {
  if (effect_type == "base") {
    variable <- "y" # shouldn't this be always one as it is the functional intercept?
  }

  if (effect_type == "binary") {
    variable <- "binary_variable"
  }

  if (effect_type == "linear") {
    variable <- "linear_variable"
  }

  if (effect_type == "smooth") {
    variable <- "smooth_variable"
  }
  unique_covariates <- unique(covariates[, variable]) # TODO check whether this works as expected
  unique_covariates
}


get_coverage_density <- function(model, theta_diff,base_range,
                         covariates, n_groups, sp, estimated_spline_densities, n_splines, X,
                         quantiles, coverage_base, coverage_binary, coverage_linear,
                         coverage_smooth) {
  alpha <- 0.05
  for (i in 1:nrow(covariates)){
    covariate_combination <- covariates[i,]
    # base
    A_base <- coverage_base$A[[1]]
    # binary
    A_binary <- coverage_binary$A[[covariates$binary_variable]]
    # linear
    A_linear <- coverage_linear$A[[covariates$linear_variable]]
    # smooth
    A_smooth <- coverage_smooth$A[[covariates$smooth_variable]]
  }
  # I don't need the range of X thetas since we only need the thetas corresponding to the y direction
  theta_diff <- theta_diff[param_range[1]:param_range[2]]
  X_for_effect <- X[, param_range[1]:param_range[2]]
  # do I need the basis with constraint or without?
  unity_matrix <- diag(1, nrow = n_splines - 1)
  # TODO check whether this produces the basis fo reach observation individually
  basis_effect <- get_basis_for_effect(effect_type, unique_covariates_for_effect) # I want the coverage for all unique covariates
  basis_functional_intercept <- X[1, base_range[0]:base_range[1]]
  S <- get_subsetting_matrix_S(param_range, length(theta_diff))

  # for the individual effects I don't need the S_j matrix directly as I can just subset V_c or V_p respectively
  # if I want to construct the KIs for the whole density I will have to construct it explicitly
  # therefore I construct it here then I can use it also in the KIs of the individual effects (not only for thetas)
  # without intercepts per covariate combination (here: 1)
  # I'm only allowed to include theta in y-direction
  if (is.null(sp)) {
    Vc <- model$Vc[(n_groups + 1 + param_range[1]):param_range[2], (n_groups + 1 + param_range[1]):param_range[2], drop = FALSE]
    Vc_inv <- try(solve(Vc), silent = TRUE)
  } else {
    Vc_inv <-  NA
  }

  Vp <- model$Vp[(n_groups + 1 + param_range[1]):param_range[2], (n_groups + 1 + param_range[1]):param_range[2], drop = FALSE]
  Vp_inv <- try(solve(Vp), silent = TRUE)
  success <- ifelse("try-error" %in% union(class(Vc_inv), class(Vp_inv)), FALSE, TRUE)

  check_coverage_Vc <- rep(NA, length(basis_effect))
  chi_statistic_Vc <- rep(NA, length(basis_effect))
  chi_statistic_Vp <- rep(NA, length(basis_effect))
  check_coverage_Vp <- rep(NA, length(basis_effect))
  A_for_effect_bases <- hashmap()
  ki_infos <- list()

  for (i in 1:length(basis_effect)) {
    mixed_basis <- kronecker(t(basis_effect[i]), unity_matrix)
    A_for_effect_base <- mixed_basis %*% S
    # I store the different As in hash table to retrieve them fast for overall density
    # KIs and coverage
    A_for_effect_bases[[unique_covariates_for_effect[i]]] <- A_for_effect_base
    Vp_mixed <- mixed_basis %*% Vp_inv %*% t(mixed_basis)
    Vp_mixed_inv <- try(solve(Vp), silent = TRUE)
    success <- ifelse("try-error" %in% union(class(Vc_inv), class(Vp_inv)), FALSE, TRUE)

    if (!any(is.na(Vc_inv))) {
      Vc_mixed <- mixed_basis %*% Vc_inv %*% t(mixed_basis)
      Vc_mixed_inv <- try(solve(Vc_mixed), silent = TRUE)
      chi_statistic_Vc[i] <- t(theta_diff) %*% Vc_mixed_inv %*% theta_diff
    } else {
      chi_statistic_Vc[i] <- NA
    }

    check_coverage_Vc[i] <- as.numeric(chi_statistic_Vc) <= qchisq(1 - alpha, df = n_splines - 1)
    chi_statistic_Vp[i] <- t(theta_diff) %*% Vp_mixed_inv %*% theta_diff
    check_coverage_Vp[i] <- as.numeric(chi_statistic_Vp) <= qchisq(1 - alpha, df = n_splines - 1)
    # calcualte KIs
    ki_infos[i] <- get_kis(estimated_spline_densities, basis_functional_intercept,  Vc_mixed, Vp_mixed, quantiles)

  }
}


get_basis_for_effect <- function(effect_type, unique_covariate_for_effect, knots) {
  for (covariate in unique_covariate_for_effect)
  if (effect_type == "base") {
    base_effect <- 1
  }

  if (effect_type == "binary") {
    base_effect <- unique_covariate_for_effect
  }

  if (effect_type == "linear") {
    base_effect <- unique_covariate_for_effect
  }

  if (effect_type == "smooth") {
    base_effect <- sum_constrained_spline_design_matrix(unique_covariate_for_effect, knots)
  }
  return(base_effect)
}


get_subsetting_matrix_S <- function(param_range, n_params) {
  dimension_unity_matrix <- param_range[2] - param_range[1]
  middle_unity_matrix <- diag(1, dimension_unity_matrix, dimension_unity_matrix)
  if (param_range[1] == 1) {
    right_zero_matrix <- matrix(0, dimension_unity_matrix, n_params - param_range[2])
    S <- cbind(middle_unity_matrix, right_zero_matrix)
  } else if (param_range[1] == n_params) {
    left_zero_matrix <- matrix(0, dimension_unity_matrix, param_range[1] - 1)
    S <- cbind(left_zero_matrix, middle_unity_matrix)
  } else {
    left_zero_matrix <- matrix(0, dimension_unity_matrix, param_range[1] - 1)
    right_zero_matrix <- matrix(0, dimension_unity_matrix, n_params - param_range[2])
    S <- cbind(left_zero_matrix, middle_unity_matrix, right_zero_matrix)
  }
  return(S)
}


get_kis <- function(spline_densities, estimated_spline_densities,
                    basis_functional_intercept,  Vc_mixed, Vp_mixed, effect_type,
                    quantiles) {
  f_hat_clr <- estimated_spline_densities$clr_density_function(quantiles, effect_type)
  f_hat <- estimated_spline_densities$density_function(quantiles)

  se_p <- sqrt(t(basis_functional_intercept) %*% Vp_mixed %*% basis_functional_intercept) # what X should I take here?
  if (!is.null(Vc_mixed)) {
    se_c <- sqrt(sqrt(t(basis_functional_intercept) %*% Vc_mixed %*% basis_functional_intercept))
  } else {
    se_c <- NA
  }
  CI_up_p <- f_hat_clr + qnorm(1 - alpha / 2) * se_p
  CI_low_p <- f_hat_clr - qnorm(1 - alpha / 2) * se_p
  CI_up_c <- f_hat_clr + qnorm(1 - alpha / 2) * se_c
  CI_low_c <- f_hat_clr - qnorm(1 - alpha / 2) * se_c
  f_true_clr <- spline_densities$clr_density_function(quantiles)
  CI_p_check <- (CI_low_p <= f_true_clr) & (f_true_clr <= CI_up_p)
  CI_c_check <- (CI_low_c <= f_true_clr) & (f_true_clr <= CI_up_c)
  CIs <- data.frame(f_true_clr, CI_low_c, CI_up_c, CI_c_check, CI_low_p, CI_up_p, CI_p_check)
  CIs
}


get_density_data <- function(densities, unpenalized, knots,
                             step_size, n_samples, sample_mode,
                             n_splines, order) {
  grid_hist <- seq(from = 0, to = 1, by = step_size)
  mids_dens <- grid_hist[1:(length(grid_hist) - 1)] + step_size / 2

  if (unpenalized) {
    support_check <- FALSE
    while (!support_check) {
      counts_dens <- sample_from_density(densities = densities,
                          n_samples = n_samples,
                          bins = grid_hist,
                          sample_mode = sample_mode,
                          quantiles = mids_dens)

      support_check <- check_knot_support(n_splines,
                                          counts_dens,
                                          mids_dens,
                                          knots,
                                          order,
                                          3)
    }
  } else {
    counts_dens <- sample_from_density(densities = densities,
                                       n_samples = n_samples,
                                       bins = grid_hist,
                                       sample_mode = sample_mode,
                                       quantiles = mids_dens)
  }


  dta_dens <- as.data.frame(cbind(counts_dens, mids_dens))
  colnames(dta_dens) <- c("counts", "y")

  Delta <- rep(step_size, length(mids_dens))

  return(list(df = dta_dens, Delta = Delta))
}

get_density_data_with_covariates <- function(densities, unpenalized, knots,
                                             step_size, sample_mode,
                                             covariates,
                                             order) {
  grid_hist <- seq(from = 0, to = 1, by = step_size)
  mids_dens <- grid_hist[1:(length(grid_hist) - 1)] + step_size / 2
  # get data frame with unique cov combis and number of replicates.
  covariates <- data.table(covariates)
  unique_covariates <- covariates[,.N, by = names(covariates)]
  n_unique_cov_combis <- nrow(unique_covariates)
  dta_dens <- data.frame(matrix(ncol = 4 + ncol(unique_covariates[,.SD, .SDcols = !c("N")]),
                                nrow = 0))
  colnames(dta_dens) <- c("counts", "y", colnames(unique_covariates[,.SD, .SDcols = !c("N")]), "Delta", "group_id")

  for (i in 1:n_unique_cov_combis) {
    n_samples <- unique_covariates[i, N]
    density <- densities[[i]]
    counts_dens <- sample_from_density(densities = density,
                                       n_samples = n_samples,
                                       bins = grid_hist,
                                       sample_mode = sample_mode,
                                       quantiles = mids_dens)
    covariates_for_observation <-  unique_covariates[replicate(length(mids_dens), i),.SD, .SDcols = !c("N")]
    # here I assume again that Delta is the same for each density component
    Delta <- rep(step_size, length(mids_dens))
    dta_dens_obs <- as.data.frame(cbind(counts_dens, mids_dens, covariates_for_observation, Delta, replicate(length(mids_dens), i)))
    colnames(dta_dens_obs) <- c("counts", "y", colnames(unique_covariates[,.SD, .SDcols = !c("N")]), "Delta", "group_id")
    dta_dens <- rbind(dta_dens, dta_dens_obs)

  }


  return(list(df = dta_dens, Delta = Delta, n_unique_cov_combis = n_unique_cov_combis, quantiles = mids_dens))
}


calculate_mse <- function(spline_densities, diff_spline_densities, norm_true) {
  # right now only the mse for the whole density can be calculated not for each effect
  relMSE <- list()
  MSE <- list()
  for (i in 1:length(spline_densities)) {
    MSE[i] <- diff_spline_densities[[i]]$norm_clr_density
    if (is.null(norm_true)) {
      norm_true <- spline_densities[[i]]$norm_clr_density
    }
    relMSE[i] <- MSE / norm_true
  }
  list("relMSE" = relMSE, "MSE" = MSE)
}


check_knot_support <- function(n_splines,
                               counts_dens,
                               mids_dens,
                               knots,
                               order,
                               threshold=0) {
  # it seems that even small counts are problematic with big N
  # therfore I introduce a threshold to mitigte this behavior
  check <- sapply(seq_len(n_splines),
                  function(k) sum(counts_dens[which(mids_dens >= knots[k] &
                                                      mids_dens <= knots[k + order])]))
  all(check > threshold)
}


plot_interpolated_density <- function(grid, # and here
                                      interpolated_density_results,
                                      true_density,
                                      scenario_number,
                                      a,
                                      b,
                                      interpolation_values,
                                      is_clr = FALSE) {
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


sample_covariates <- function(n_obs) {
  binary <- sample(c(0,1), n_obs, replace = TRUE)
  linear <- rnorm(n_obs) * 100
  smooth <- rnorm(n_obs) * 5 + 5
  covariates <- data.frame(
    binary_variable = binary,
    linear_variable = linear,
    smooth_variable = smooth
  )

  covariates
}


smooth.construct.md.smooth.spec <- function(object, data, knots) {
  x <- data[[object$term]]
  # getting specifications for continuous component
  if (length(object$p.order) == 1)
    m <- rep(object$p.order, 2)
  else m <- object$p.order  # m[1] - basis order, m[2] - penalty order
  m[is.na(m)] <- 2 ## default
  object$p.order <- m
  if (object$bs.dim < 0)
    object$bs.dim <- max(10, m[1] + 1) ## default
  nk <- object$bs.dim - m[1] # basis dimension - order of spline -> number of interior knots for continuous component
  if (nk <= 0)
    stop("Basis dimension too small for b-spline order")
  if (length(object$term) != 1)
    stop("Basis only handles 1D smooths")
  cont_dim <- object$bs.dim
  object$xt <- object$xt[[1]]
  # initializing discrete component
  if (is.null(object$xt$values_discrete)) {
    object$xt$values_discrete <- c(0, 1)
  }
  t_discrete <- object$xt$values_discrete
  if (is.null(object$xt$weights_discrete)) {
    object$xt$weights_discrete <- rep(1, length(t_discrete))
  }
  w_discrete <- object$xt$weights_discrete
  if (length(t_discrete) != length(w_discrete)) {
    stop("Lengths of values_discrete and weights_discrete have to be the same.")
  }

  if (is.null(object$xt$domain_continuous)) {
    object$xt$domain_continuous <- range(t_discrete, x)
  }
  cont_positions <- which(!(x %in% t_discrete))
  x_cont <- x[cont_positions]
  if (object$xt$domain_continuous[1] > min(x_cont) ||
      object$xt$domain_continuous[2] < max(x_cont)) {
    stop("Given domain does not include data corresponding to continuous component.")
  }

  # add discrete value t_{D+1} and weight corresponding to the continuous component
  t_discrete[length(t_discrete) + 1] <- range(t_discrete, x)[2] + 1
  w_discrete[length(w_discrete) + 1] <- diff(object$xt$domain_continuous)

  discrete_dim <- length(t_discrete)
  object$bs.dim <- object$bs.dim + discrete_dim # combined dimension (before implementing constraints!)
  # set all continuous values to t_{D+1} for discrete basis
  x_discrete <- x
  x_discrete[cont_positions] <- max(t_discrete)

  k <- sort(knots[[object$term]])
  if (is.null(k)) {
    xl <- object$xt$domain_continuous[1]
    xu <- object$xt$domain_continuous[2]
  } else if (length(k) == 2) {
    xl <- min(k)
    xu <- max(k)
    if (xl > min(x_cont) || xu < max(x_cont))
      stop("Knot range does not include data corresponding to continuous component.")
    # if (xl == 0 || xu == 1) {
    #   stop("Knots 0 and 1 are reserved for discrete
    #        spline part and shouldn't be provided")
    # } ### E: No, see below
  }
  if (is.null(k) || length(k) == 2) {
    # note that the xl > 0 and xu < 1 at all times and the same must
    # hold for all knots, as 0 and 1 are knots for the discrete part
    ### E: No, the knots can also be <= 0 or >= 1, but the resulting splines are
    ###    only evaluated on (0, 1), here x_cont. Since equidistant knots are
    ###    preferable, I changed it to the knot extension used in the smooth
    ###    construct for "ps" (without shifting).
    # nk = basis dimension - order of spline
    xr <- xu - xl
    # xl <- xl - xr * 0.001 ### E: this shifts the knots slightly outwards, i.e., not (exactly) the supplied knots (for length(k) == 2) are used
    # xu <- xu + xr * 0.001 ### E: I commented it out, since for me it's unclear why not to used the supplied knots
    dx <- (xu - xl)/(nk - 1)
    k <- seq(xl - dx * (m[1] + 1), xu + dx * (m[1] + 1),
             length = nk + 2 * m[1] + 2)
  } else {
    if (length(k) != nk + 2 * m[1] + 2)
      stop(paste("There should be ", nk + 2 * m[1] + 2, " supplied knots"))
    # if (0 %in% k || 1  %in% k) {
    #   stop("Knots 0 and 1 are reserved for discrete
    #        spline part and shouldn't be provided")
    # } ### E: No, see above
  }
  ord <- m[1] + 2
  if (k[ord] != object$xt$domain_continuous[1] ||
      k[length(k) - (ord - 1)] != object$xt$domain_continuous[2])
    warning("Knots do not match domain of continuous component.")
  if (is.null(object$deriv)) {
    object$deriv <- 0
  } else if (object$deriv != 0) {
    warning("The mixed density smoother is not intended to be used for derivatives. Reasonable behavior is only guaranteed for deriv = 0.")
  }
  # construct design matrices from transformed B-splines integrating to zero (see
  # Appendix B of Maier et al., 2021, based on Wood, 2017, Section 1.8.1) and
  # apply embedding to combine them to one design matrix of a mixed basis (Maier
  # et al., 2021, Proposition A.4; Maier et al., 2022, Section 2.2)
  C_cont <- sapply(1:(length(k) - ord), function(j) integrate(function(x)
    splines::splineDesign(knots = k, x, ord = ord, derivs = object$deriv)[, j],
    lower = k[ord], upper = k[length(k) - (ord - 1)])$value)
  Z_cont <- MASS::Null(C_cont)
  design_cont <- matrix(0, nrow = length(x), ncol = ncol(Z_cont))
  design_cont[cont_positions, ] <- splines::splineDesign(knots = k, x_cont, ord = ord,
                                                         derivs = object$deriv) %*% Z_cont

  k_discrete <- sapply(seq_len(length(t_discrete) + 1),
                       function(j) mean(c(min(t_discrete) - 1, t_discrete, max(t_discrete) + 1)[j:(j+1)]))
  object$knots_discrete <- k_discrete
  object$values_discrete <- t_discrete
  C_discrete <- object$weights_discrete <- w_discrete # integrals of basis functions are equal to weights
  Z_discrete <- MASS::Null(C_discrete)
  design_discrete <- splines::splineDesign(k_discrete, x_discrete, 1)  %*% Z_discrete

  object$X <- cbind(design_cont, design_discrete) # combine both design matrices to get final design matrix
  object$Z <- list(cont = Z_cont, discrete = Z_discrete)

  if (!is.null(k)) {
    if (sum(colSums(object$X) == 0) > 0)
      warning("There is *no* information about some basis coefficients")
  }
  if (length(unique(x)) < object$bs.dim)
    warning("Basis dimension is larger than number of unique covariates")
  if (is.null(object$mono))
    object$mono <- 0
  if (object$mono != 0) {
    stop("SCOP splines are not supported yet!")
  } else {
    # construct penalty matrices including necessary transformation (see Appendix
    # B of Maier et al., 2021) and combine them to one penalty (block) matrix
    # for the mixed design matrix
    if (m[2] > 0) {
      D_cont <- diff(diag(cont_dim), differences = m[2])
    } else {
      D_cont <- diag(cont_dim)
    }
    S_cont <- t(Z_cont) %*% crossprod(D_cont) %*% Z_cont

    if (is.null(object$xt$penalty_discrete)) {
      D_discrete <- matrix(0, nrow = discrete_dim, ncol = discrete_dim) # default: discrete component unpenalized
    } else if (object$xt$penalty_discrete > 0) {
      D_discrete <- diff(diag(discrete_dim), differences = object$xt$penalty_discrete)
    } else if (object$xt$penalty_discrete == 0) {
      D_discrete <- diag(discrete_dim)
    } else {
      warning("penalty_discrete has to be a non-negative integer. No penalty for discrete component is used.")
      D_discrete <- matrix(0, nrow = discrete_dim, ncol = discrete_dim)
    }
    S_discrete <- t(Z_discrete) %*% crossprod(D_discrete) %*% Z_discrete

    # object$D <- list(continuous = D_cont, discrete = D_discrete)
    # I decided to comment out the definition of object$D since it is not even
    # listed under Value in ?smooth.construct

    # combine penalties in a block matrix
    S <- cbind(rbind(S_cont, matrix(0, ncol = ncol(S_cont), nrow = nrow(S_discrete))),
               rbind(matrix(0, nrow = nrow(S_cont), ncol = ncol(S_discrete)), S_discrete))

    object$S <- list(S)
    object$rank <- ifelse(is.null(object$xt$penalty_discrete), cont_dim - m[2],
                          ifelse(object$xt$penalty_discrete == 0,
                                 cont_dim - m[2] + discrete_dim - 1, # we loose one dimension when applying constraint
                                 cont_dim - m[2] + discrete_dim - object$xt$penalty_discrete))
    object$null.space.dim <- ifelse(is.null(object$xt$penalty_discrete), m[2],
                                    m[2] + object$xt$penalty_discrete)
  }
  object$knots <- k
  object$m <- m
  class(object) <- "mdspline.smooth"
  object
}


smooth.construct.d.smooth.spec <- function(object, data, knots) {
  x <- data[[object$term]]
  object$xt <- object$xt[[1]]
  # getting specifications for continuous component
  if (length(object$p.order) == 1)
    m <- rep(object$p.order, 2)
  else m <- object$p.order  # m[1] - basis order, m[2] - penalty order
  m[is.na(m)] <- 2 ## default
  object$p.order <- m
  if (object$bs.dim < 0)
    object$bs.dim <- max(10, m[1] + 1) ## default
  nk <- object$bs.dim - m[1] # basis dimension - order of spline -> number of interior knots for continuous component
  if (nk <= 0)
    stop("Basis dimension too small for b-spline order")
  if (length(object$term) != 1)
    stop("Basis only handles 1D smooths")
  cont_dim <- object$bs.dim
  if (is.null(object$xt$domain_continuous)) {
    object$xt$domain_continuous <- c(0, 1)
  }
  x_cont <- x
  if (object$xt$domain_continuous[1] > min(x_cont) ||
      object$xt$domain_continuous[2] < max(x_cont)) {
    stop("Given domain does not include data corresponding to continuous component.")
  }

  k <- sort(knots[[object$term]])
  if (is.null(k)) {
    xl <- object$xt$domain_continuous[1]
    xu <- object$xt$domain_continuous[2]
  } else if (length(k) == 2) {
    xl <- min(k)
    xu <- max(k)
    if (xl > min(x_cont) || xu < max(x_cont))
      stop("Knot range does not include data corresponding to continuous component.")
    # if (xl == 0 || xu == 1) {
    #   stop("Knots 0 and 1 are reserved for discrete
    #        spline part and shouldn't be provided")
    # } ### E: No, see below
  }
  if (is.null(k) || length(k) == 2) {
    # note that the xl > 0 and xu < 1 at all times and the same must
    # hold for all knots, as 0 and 1 are knots for the discrete part
    ### E: No, the knots can also be <= 0 or >= 1, but the resulting splines are
    ###    only evaluated on (0, 1), here x_cont. Since equidistant knots are
    ###    preferable, I changed it to the knot extension used in the smooth
    ###    construct for "ps" (without shifting).
    # nk = basis dimension - order of spline
    xr <- xu - xl
    # xl <- xl - xr * 0.001 ### E: this shifts the knots slightly outwards, i.e., not (exactly) the supplied knots (for length(k) == 2) are used
    # xu <- xu + xr * 0.001 ### E: I commented it out, since for me it's unclear why not to used the supplied knots
    dx <- (xu - xl)/(nk - 1)
    k <- seq(xl - dx * (m[1] + 1), xu + dx * (m[1] + 1),
             length = nk + 2 * m[1] + 2)
  } else {
    if (length(k) != nk + 2 * m[1] + 2)
      stop(paste("There should be ", nk + 2 * m[1] + 2, " supplied knots"))
    # if (0 %in% k || 1  %in% k) {
    #   stop("Knots 0 and 1 are reserved for discrete
    #        spline part and shouldn't be provided")
    # } ### E: No, see above
  }
  ord <- m[1] + 2
  if (k[ord] != object$xt$domain_continuous[1] ||
      k[length(k) - (ord - 1)] != object$xt$domain_continuous[2])
    warning("Knots do not match domain of continuous component.")
  if (is.null(object$deriv)) {
    object$deriv <- 0
  } else if (object$deriv != 0) {
    warning("The density smoother is not intended to be used for derivatives. Reasonable behavior is only guaranteed for deriv = 0.")
  }
  # construct design matrices from transformed B-splines integrating to zero (see
  # Appendix B of Maier et al., 2021, based on Wood, 2017, Section 1.8.1) and
  # apply embedding to combine them to one design matrix of a mixed basis (Maier
  # et al., 2021, Proposition A.4; Maier et al., 2022, Section 2.2)
  C_cont <- sapply(1:(length(k) - ord), function(j) integrate(function(x)
    splines::splineDesign(knots = k, x, ord = ord, derivs = object$deriv)[, j],
    lower = k[ord], upper = k[length(k) - (ord - 1)])$value)
  Z_cont <- MASS::Null(C_cont)
  design_cont <- splines::splineDesign(knots = k, x_cont, ord = ord,
                                       derivs = object$deriv) %*% Z_cont

  object$X <- design_cont # combine both design matrices to get final design matrix
  object$Z <- list(cont = Z_cont)

  if (!is.null(k)) {
    if (sum(colSums(object$X) == 0) > 0)
      warning("There is *no* information about some basis coefficients")
  }
  if (length(unique(x)) < object$bs.dim) # evtl -1, wegen integrate to zero constraint
    warning("Basis dimension is larger than number of unique covariates")
  if (is.null(object$mono))
    object$mono <- 0
  if (object$mono != 0) {
    stop("SCOP splines are not supported yet!")
  } else {
    # construct penalty matrices including necessary transformation (see Appendix
    # B of Maier et al., 2021) and combine them to one penalty (block) matrix
    # for the mixed design matrix
    if (m[2] > 0) {
      D_cont <- diff(diag(cont_dim), differences = m[2])
    } else {
      D_cont <- diag(cont_dim)
    }
    S_cont <- t(Z_cont) %*% crossprod(D_cont) %*% Z_cont

    # object$D <- list(continuous = D_cont, discrete = D_discrete)
    # I decided to comment out the definition of object$D since it is not even
    # listed under Value in ?smooth.construct

    # combine penalties in a block matrix
    S <- S_cont

    object$S <- list(S)
    object$rank <- cont_dim - m[2] - 1
    object$null.space.dim <- m[2]
  }
  object$knots <- k
  object$m <- m
  class(object) <- "dspline.smooth"
  object
}


Predict.matrix.dspline.smooth <- function(object, data) {
  m <- object$m[1] + 1
  ll <- object$xt$domain_continuous[1] # object$knots[m + 1]
  ul <- object$xt$domain_continuous[2] # object$knots[length(object$knots) - m]
  m <- m + 1
  x <- data[[object$term]]
  x_cont <- x
  ind <- list(cont = (x_cont <= ul & x_cont >= ll))
  if (is.null(object$deriv)) {
    object$deriv <- 0
  } else if (object$deriv != 0) {
    warning("The mixed density smoother is not intended to be used for derivatives. Reasonable behavior is only guaranteed for deriv = 0.")
  }
  if (sum(sapply(ind, sum)) == length(x)) {
    k <- object$knots
    Z_cont <- object$Z$cont
    design_cont <- splines::splineDesign(knots = k, x_cont, ord = m,
                                                           derivs = object$deriv) %*% Z_cont
    X <- design_cont
  } else {
    stop("Supplied data is not in support of underlying density!")
  }
  if (object$mono == 0){
    X
  } else {
    stop("SCOP splines are not supported yet!")
  }
}





