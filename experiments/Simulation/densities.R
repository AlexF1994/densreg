library("purrr")
library("truncnorm")

eps = 0.000001


get_densities <- function(density_params, calculate_norm=FALSE, param_scale=0,
                          component=FALSE, ...) {
  implemented_densities <- c("constant", "line", "beta", "truncated_normal", "spline")

  if (density_params$density_name == "constant") {
    density_function <- partial(constant_density, constant = 1)
    distribution_function <- NULL
  }

  if (density_params$density_name == "line") {
    density_function <- partial(line_density, slope = density_params$slope)
    distribution_function <- partial(line_distribution, slope = density_params$slope)
  }

  if (density_params$density_name == "beta") {
    # take absolute value to avoid negative parmater values for a and b
    a <- ifelse(param_scale, (1 + abs(param_scale)) * density_params$a, density_params$a)
    b <- ifelse(param_scale, (1 + abs(param_scale)) * density_params$b, density_params$b)
    density_function <- partial(beta_density, a = a,
                                b = b)
    distribution_function <- partial(beta_distribution, a = a,
                                b = b)
  }

  if (density_params$density_name == "truncated_normal") {
    mean <- ifelse(param_scale, param_scale * density_params$mean, density_params$mean)
    sd <- ifelse(param_scale, param_scale * density_params$sd, density_params$sd)
    density_function <- partial(truncated_normal_density,
                                mean = mean,
                                sd = sd)
    distribution_function <- partial(truncated_normal_distribution,
                                mean = mean,
                                sd = sd)
  }

  if (density_params$density_name == "spline") {
    # clr_density has to be declared seperately for spline densities
    # as otherwise we have an integration in integration problem
    # which takes ages to compute
    # I assume that the relevant theta parts for every density component
    # are always at the same position of the theta matrix
    if (!is.null(component)) {
      theta <- get_theta_for_component(density_params, component)
    } else {
      theta <- density_params$theta
    }

    if (component == "smooth") {
      kwargs <- list(...)
      clr_density <- get_multivariate_spline_clr_density(theta,
                                            density_params$knots, # for now I assume the same knots for every component
                                            density_params$knots_smooth_covariate,
                                            density_params$order,
                                            kwargs$smooth_variable)
    } else{
    clr_density <- get_spline_clr_density(theta,
                                          density_params$knots, # for now I assume the same knots for every component
                                          density_params$order)
    }

    density_function <- partial(spline_density,
                                spline_clr_density = clr_density)
    distribution_function <- NULL
  }

  if (density_params$density_name != "spline") {
    if (density_params$density_name == "constant") {
      clr_density <- partial(constant_density, constant = 0)
    } else {
      clr_density <- partial(clr, density_function = density_function)
    }
  }

  if (calculate_norm) {
    norm_clr_density <- integrate(partial(clr_density_squared,
                                          clr_density = clr_density),
                                  lower = 0 + eps,
                                  upper = 1 - eps, subdivisions = 1000)$value
  }
  else {
    norm_clr_density <- NULL
  }

  return(
    list(density_function = density_function,
         clr_density_function = clr_density,
         norm_clr_density = norm_clr_density,
         density_name = density_params$density_name,
         distribution_function = distribution_function
         )
    )
}


get_theta_for_component <- function(density_params, component) {
  # theta is assumed to be without intercept(s)
  complete_theta <- density_params$theta
  if (component == "base") {
    theta_range <- density_params$base_range
  }

  if (component == "binary") {
    theta_range <- density_params$binary_range
  }

  if (component == "linear") {
    theta_range <- density_params$linear_range
  }

  if (component == "smooth") {
    theta_range <- density_params$smooth_range
  }
  theta_for_component <- complete_theta[theta_range[1]:theta_range[-1]]
  return(theta_for_component)
}


get_densities_with_covariates <- function(density_params, covariates, calculate_norm=FALSE) {
  # I get one density for every observation
  # I will order the densities via observation number --> list with n entries
  n_obs <- nrow(covariates)
  densities <- lapply(1:n_obs, get_densities_with_covariates_single,
                      covariates = covariates, density_params = density_params,
                      calculate_norm = calculate_norm)

  return(densities)
}

get_densities_with_covariates_single <- function(observation_index, density_params, covariates, calculate_norm) {
  # base density
  base_density_component <- get_densities(density_params, calculate_norm = FALSE, component = "base")
  # binary density
  if (covariates$binary_variable[observation_index] == 1) {
    binary_density_component <- get_densities(density_params, calculate_norm = FALSE, component = "binary")
  }

  else {
    binary_density_component <- get_densities(list(density_name = "constant"), calculate_norm = FALSE)
  }

  # linear density
  linear_variable <- covariates$linear_variable[observation_index]
  density_part <- get_densities(density_params, calculate_norm = FALSE, component = "linear")
  linear_density_component <- list(density_function = partial(pertubated_density, density_function = density_part$density_function,
                                                              pertubator = linear_variable),
                                   clr_density_function =  partial(scaled_clr_density, clr_density_function = density_part$clr_density_function,
                                                                   pertubator = linear_variable))
  # smooth density
  smooth_variable = covariates$smooth_variable[observation_index]


  smooth_density_component <- get_densities(density_params, param_scale = smooth_variable, component = "smooth", smooth_variable = smooth_variable)

  density_components <- list(base_density_component,
                             binary_density_component,
                             linear_density_component,
                             smooth_density_component)

  density_function <- partial(composite_density, components = density_components)
  clr_density_function <- partial(composite_clr_density, components = density_components)

  densities <- list(density_function = density_function,
                    clr_density_function = clr_density_function,
                    density_name = density_params$density_name
  )

  if (calculate_norm) {
    norm_clr_density <- integrate(partial(clr_density_squared,
                                          clr_density = clr_density_function),
                                  lower = 0 + eps,
                                  upper = 1 - eps,
                                  subdivisions = 1000)$value
    densities <- list.append(densities, norm_clr_density = norm_clr_density)
  }

  return(
    densities
  )

}


# true density to interpolate (we use a beta distribution)
beta_density <- function(quantiles, a, b) {
  dbeta(quantiles, a, b)
}

beta_distribution <- function(quantiles, a, b) {
  pbeta(quantiles, a, b)
}


line_density <- function(quantiles, slope = 0.1) {
  b <- 1 - 0.5 * slope
  return(b + slope * quantiles)
}

constant_density <- function(quantiles, constant = 1) {
  return(replicate(length(quantiles), constant))
}

line_distribution <- function(quantiles, slope = 0.1) {
  b <- 1 - 0.5 * slope
  return(b * quantiles + (slope / 2) * quantiles^2)
}


truncated_normal_density <- function(quantiles, mean=0, sd=1) {
  dtruncnorm(quantiles, 0, 1, mean, sd)
}

truncated_normal_distribution <- function(quantiles, mean=0, sd=1) {
  ptruncnorm(quantiles, mean, sd)
}

pertubated_density <- function(density_function, pertubator, quantiles) {
  return(density_function(quantiles)^pertubator)
}

scaled_clr_density <- function(clr_density_function, pertubator, quantiles) {
  return(clr_density_function(quantiles) * pertubator)
}

composite_density <- function(components, quantiles, component_name = "all") {
  if (component_name == "all") {
    density_values <- c(rep(1, length(quantiles)))
    for (component in components) {
      density_values <- density_values * component$density_function(quantiles)
    }
    density_values
  } else if (component_name == "base") {
    components[[1]]$density_function(quantiles)
  } else if (component_name == "binary") {
    components[[2]]$density_function(quantiles)
  } else if (component_name == "linear") {
    components[[3]]$density_function(quantiles)
  } else if (component_name == "smooth") {
    components[[4]]$density_function(quantiles)
  }
}

composite_clr_density <- function(components, quantiles, component_name = "all") {
  if (component_name == "all") {
    clr_density_values <- c(rep(0, length(quantiles)))
    for (component in components) {
      clr_density_values <- clr_density_values + component$clr_density_function(quantiles)
    }
    clr_density_values
  } else if (component_name == "base") {
    components[[1]]$clr_density_function(quantiles)
  } else if (component_name == "binary") {
    components[[2]]$clr_density_function(quantiles)
  } else if (component_name == "linear") {
    components[[3]]$clr_density_function(quantiles)
  } else if (component_name == "smooth") {
    components[[4]]$clr_density_function(quantiles)
  }
}

clr <- function(density_function, quantiles) {
  integral_log_density <- integrate(function(x) log(density_function(x)),
                           lower = 0 + eps, upper = 1 - eps, subdivisions = 1000)$value
  log(density_function(quantiles)) - integral_log_density
}


clr_density_squared <- function(clr_density, quantiles) {
  clr_density(quantiles)^2
}


inverse_clr <- function(clr_density, quantiles){
  integral_exp_clr <- integrate(function(x) exp(clr_density(x)),
                                lower = 0 + eps, upper = 1 - eps, subdivisions = 1000)$value

  exp(clr_density(quantiles)) / integral_exp_clr
}

get_spline_clr_density <- function(theta, knots, order) {
  partial(spline_clr_density, theta = theta, knots = knots, order = order)
}


get_multivariate_spline_clr_density <- function(theta, knots, knots_covariate, order, covariate) {
  partial(multivariate_spline_clr_density, theta = theta, knots = knots,
          knots_covariate = knots_covariate, order = order, covariate = covariate)
}


spline_clr_density <- function(quantiles, theta, knots, order) {

  design_matrix <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = order)
  values <- design_matrix %*% theta

  values
}


multivariate_spline_clr_density <- function(quantiles, theta, knots, knots_covariate, order, covariate) {

  design_matrix_y <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = order)
  design_matrix_x <- sum_constrained_spline_design_matrix(x = covariate, knots = knots_covariate, ord = order)
  design_matrix <- kronecker(design_matrix_x, design_matrix_y)
  values <- design_matrix %*% theta

  values
}


spline_density <- function(quantiles, spline_clr_density) {
  inverse_clr(spline_clr_density, quantiles = quantiles)
}


# constrained_spline_design_matrix creates a B-spline basis transformed to fulfill the integrate-to-zero
# constraint. Arguments are the same as in splines:splineDesign:
# - x: a numeric vector of values at which to evaluate the basis functions or
#   derivatives. The values in x must be between the “inner” knots knots[ord] and
#   knots[ length(knots) - (ord-1)].
# - knots: a numeric vector of (inner and outer) knot positions (which will be
#   sorted increasingly if needed)
# - ord: a positive integer giving the order of the spline function. This is the
#   number of coefficients in each piecewise polynomial segment, thus a cubic
#   spline has order 4. Defaults to 4.
# Value: A matrix with length(x) rows and length(knots) - ord - 1 columns. The i'th
# row of the matrix contains the evaluations of the constrained spline functions
# (defined by the knot vector and the order) at the i'th value of x.

constrained_spline_design_matrix <- function(x, knots, ord = 4) {
  C <- sapply(1:(length(knots) - ord), function(j) integrate(function(x)
    splines::splineDesign(knots = knots, x, ord = ord, derivs = 0)[, j],
    lower = knots[ord], upper = knots[length(knots) - (ord-1)])$value)
  Z <- MASS::Null(C)
  X_L20 <- splines::splineDesign(knots = knots, x, ord = ord, derivs = 0) %*% Z
}

constrained_tensor_spline_design_matrix <- function(x, knots, ord = 4) {
  C <- sapply(1:(length(knots) - ord), function(j) integrate(function(x)
    splines::splineDesign(knots = knots, x, ord = ord, derivs = 0)[, j],
    lower = knots[ord], upper = knots[length(knots) - (ord-1)])$value)
  Z <- MASS::Null(C)
  X_L20 <- splines::splineDesign(knots = knots, x, ord = ord, derivs = 0) %*% Z
}

sum_constrained_spline_design_matrix <- function(x, knots, ord = 4) {
  X <- splines::splineDesign(knots = knots, x)
  C <- rep(1, nrow(X)) %*% X
  qrc <- qr(t(C))
  Z <- qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]
  X_0 <- X %*% Z
  X_0
}


sample_from_density <- function(densities, n_samples, bins,
                                sample_mode = c("bin", "value"),
                                quantiles = NULL) {
  sample_mode <- match.arg(sample_mode)
  n_bins <- length(bins) - 1
  probs <- get_bin_probabilities(bins, densities, quantiles)

  if (sample_mode == "value") {
    samples <- rep(NA, n_samples)
    for (i in seq_len(n_samples)) {
      random_bin <- sample(n_bins, 1, prob = probs)
      samples[i] <- runif(1, min = bins[random_bin], max = bins[random_bin + 1])
    }
  }
  else {
    samples <- rmultinom(n = 1, size = n_samples, prob = probs)
  }

  samples

}

# should get rid of integrate
# ideas:
# - use sum for spline approximations
# use quantiles for true densities
get_bin_probabilities <- function(bins, densities, quantiles=NULL) {
  eps <- 0.0000001
  n_bins <- length(bins) - 1
  if (densities$density_name == "spline") {
    step_size = 1 / n_bins
    linear_predictor <- c(log(step_size) + densities$clr_density_function(quantiles))
    bin_probabilities <- (exp(linear_predictor) + (eps / n_bins)) / (sum(exp(linear_predictor)) + eps)
    if (any(is.nan(bin_probabilities))) {
      print("here")
    }
  }
  # here I can use the quantile functions if it still takes too long to compute
  else {
    n_bins <- length(bins) - 1
    bin_probabilities <- c()
    for (i in 1:n_bins) {
      bin_probabilities[i] <- densities$distribution_function(bins[i+1]) - densities$distribution_function(bins[i])
    }
    bin_probabilities
  }
  bin_probabilities
}

