library("purrr")
library("truncnorm")
library("rlist")



get_densities <- function(density_params, calculate_norm=FALSE, param_scale=0,
                          component=FALSE, Z_cont=NULL, ...) {
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
    # take absolute value to avoid negative paramater values for a and b
    a <- ifelse(param_scale, (abs(param_scale)) + density_params$a, density_params$a)
    b <- density_params$b
    density_function <- partial(beta_density, a = a,
                                b = b)
    distribution_function <- partial(beta_distribution, a = a,
                                b = b)
  }

  if (density_params$density_name == "truncated_normal") {
    mean <- ifelse(param_scale, param_scale * density_params$mean, density_params$mean)
    sd <- density_params$sd
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
                                            kwargs$smooth_variable,
                                            Z_cont,
                                            kwargs$smooth_covariate_design_matrix)
    } else{
    clr_density <- get_spline_clr_density(theta,
                                          density_params$knots, # for now I assume the same knots for every component
                                          density_params$order,
                                          Z_cont)
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
                                  lower = 0 ,
                                  upper = 1, subdivisions = 100)$value
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


get_densities_with_covariates <- function(density_params, covariates, smooth_covariate_design_matrix=NULL,
                                          calculate_norm=FALSE) {
  # I get one density for every observation
  # I will order the densities via observation number --> list with n entries
  Z_cont <- NULL
  # this part could be also done on the highest level
  if (density_params$density_name == "spline") {
    C_cont <- sapply(1:(length(density_params$knots) - density_params$order), function(j) integrate(function(x)
      splines::splineDesign(knots = density_params$knots, x, ord = density_params$order, derivs = 0)[, j],
      lower = density_params$knots[density_params$order], upper = density_params$knots[length(density_params$knots) - (density_params$order - 1)])$value)
    Z_cont <- MASS::Null(C_cont)
  }
  n_obs <- nrow(covariates)
  densities <- lapply(1:n_obs, get_densities_with_covariates_single,
                      covariates = covariates, density_params = density_params,
                      calculate_norm = calculate_norm, Z_cont = Z_cont,
                      smooth_covariate_design_matrix = smooth_covariate_design_matrix)

  return(densities)
}

get_densities_with_covariates_single <- function(observation_index, density_params, covariates, calculate_norm, Z_cont,
                                                 smooth_covariate_design_matrix) {
  # base density
  base_density_component <- get_densities(density_params, calculate_norm = calculate_norm, component = "base", Z_cont = Z_cont)
  # binary density
  if (covariates$binary_variable[observation_index] == 1) {
    density_params_binary <- get_binary_density_params(density_params)
    binary_density_component <- get_densities(density_params_binary, calculate_norm = calculate_norm, component = "binary", Z_cont = Z_cont)
  }

  else {
    binary_density_component <- get_densities(list(density_name = "constant"), calculate_norm = calculate_norm)
  }

  # linear density
  linear_variable <- covariates$linear_variable[observation_index]
  density_params_linear <- get_linear_density_params(density_params)
  density_part <- get_densities(density_params_linear, calculate_norm = FALSE, component = "linear", Z_cont = Z_cont)
  linear_density_component <- list(density_function = partial(pertubated_density, density_function = density_part$density_function,
                                                              pertubator = linear_variable),
                                   clr_density_function =  partial(scaled_clr_density, clr_density_function = density_part$clr_density_function,
                                                                   pertubator = linear_variable))
  # smooth density
  smooth_variable = covariates$smooth_variable[observation_index]
  smooth_variable_design_matrix <- NULL

  if (!is.null(smooth_covariate_design_matrix)) {
    smooth_variable_design_matrix <- smooth_covariate_design_matrix[observation_index,]
  }


  smooth_density_component <- get_densities(density_params, param_scale = smooth_variable, component = "smooth",
                                            calculate_norm = calculate_norm,
                                            smooth_variable = smooth_variable,
                                            Z_cont = Z_cont,
                                            smooth_covariate_design_matrix = smooth_variable_design_matrix)

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
    norm_clr_density_linear <- integrate(partial(clr_density_squared,
                                          clr_density = linear_density_component$clr_density_function),
                                  lower = 0,
                                  upper = 1,
                                  subdivisions = 100)$value
    norm_clr_density <- integrate(partial(clr_density_squared,
                                          clr_density = clr_density_function),
                                  lower = 0,
                                  upper = 1,
                                  subdivisions = 100)$value
    densities <- list.append(densities,
                             norm_clr_density = norm_clr_density,
                             norm_clr_density_base = base_density_component$norm_clr_density,
                             norm_clr_density_binary = binary_density_component$norm_clr_density,
                             norm_clr_density_linear = norm_clr_density_linear,
                             norm_clr_density_smooth = smooth_density_component$norm_clr_density
                             )
  }

  return(
    densities
  )

}


get_binary_density_params <- function(density_params) {
  density_params_binary <- density_params
  if (density_params$density_name == "beta") {
    density_params_binary$b <- density_params$b + 1
  }
  else {
    density_params_binary$sd <- density_params$sd - 0.25
  }
  return(density_params_binary)
}


get_linear_density_params <- function(density_params) {
  density_params_linear <- density_params
  if (density_params$density_name == "beta") {
    density_params_linear$b <- density_params$b - 0.5
  }
  else {
    density_params_linear$sd <- density_params$sd + 0.25
  }
  return(density_params_linear)
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
                           lower = 0 , upper = 1, subdivisions = 100)$value
  log(density_function(quantiles)) - integral_log_density
}


clr_density_squared <- function(clr_density, quantiles) {
  clr_density(quantiles)^2
}


inverse_clr <- function(clr_density, quantiles){
  integral_exp_clr <- integrate(function(x) exp(clr_density(x)),
                                lower = 0 , upper = 1, subdivisions = 100)$value

  exp(clr_density(quantiles)) / integral_exp_clr
}

get_spline_clr_density <- function(theta, knots, order, Z_cont) {
  partial(spline_clr_density, theta = theta, knots = knots, order = order, Z_cont = Z_cont)
}


get_multivariate_spline_clr_density <- function(theta, knots, knots_covariate, order, covariate, Z_cont,
                                                smooth_covariate_design_matrix) {
  partial(multivariate_spline_clr_density, theta = theta, knots = knots,
          knots_covariate = knots_covariate, order = order, covariate = covariate,
          Z_cont = Z_cont, smooth_covariate_design_matrix = smooth_covariate_design_matrix)
}


spline_clr_density <- function(quantiles, theta, knots, order, Z_cont) {

  design_matrix <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = order, Z_cont = Z_cont)
  values <- design_matrix %*% theta

  values
}


multivariate_spline_clr_density <- function(quantiles, theta, knots, knots_covariate, order, covariate, Z_cont,
                                            smooth_covariate_design_matrix) {

  design_matrix_y <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = order, Z_cont = Z_cont)
  design_matrix_x <- smooth_covariate_design_matrix
  design_matrix <- kronecker(t(design_matrix_x), design_matrix_y)
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

constrained_spline_design_matrix <- function(x, knots, ord = 4, Z_cont = NULL) {
  if (is.null(Z_cont)) {
    C_cont <- sapply(1:(length(knots) - ord), function(j) integrate(function(x)
      splines::splineDesign(knots = knots, x, ord = ord, derivs = 0)[, j],
      lower = knots[ord], upper = knots[length(knots) - (ord - 1)])$value)
    Z_cont <- MASS::Null(C_cont)
  }
  design_cont <- splines::splineDesign(knots = knots, x, ord = ord,
                                                         derivs = 0) %*% Z_cont

  X_L20 <- design_cont
  X_L20
}

constrained_tensor_spline_design_matrix <- function(x, knots, ord = 4) {
  C <- sapply(1:(length(knots) - ord), function(j) integrate(function(x)
    splines::splineDesign(knots = knots, x, ord = ord, derivs = 0)[, j],
    lower = knots[ord], upper = knots[length(knots) - (ord-1)])$value)
  Z <- MASS::Null(C)
  X_L20 <- splines::splineDesign(knots = knots, x, ord = ord, derivs = 0) %*% Z
}

sum_constrained_spline_design_matrix <- function(x, knots, ord = 4) {
  X <- splines::splineDesign(knots = knots, x, ord = ord, derivs = 0)
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
  n_bins <- length(bins) - 1
  if (densities$density_name == "spline") {
    step_size = 1 / n_bins
    linear_predictor <- c(log(step_size) + densities$clr_density_function(quantiles))

    bin_probabilities <- (exp(linear_predictor) / (sum(exp(linear_predictor))))
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

