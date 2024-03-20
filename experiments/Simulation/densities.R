library("purrr")
library("truncnorm")


get_densities <- function(density_params, calculate_norm=FALSE) {
  implemented_densities <- c("line", "beta", "truncated_normal", "spline")

  if (density_params$density_name == "line") {
    density_function <- partial(line_density, slope = density_params$slope)
  }

  if (density_params$density_name == "beta") {
    density_function <- partial(beta_density, a = density_params$a,
                                b = density_params$b)
  }

  if (density_params$density_name == "truncated_normal") {
    density_function <- partial(truncated_normal_density,
                                mean = density_params$mean,
                                sd = density_params$sd)
  }

  if (density_params$density_name == "spline") {
    density_function <- partial(spline_density, density_params$theta,
                                density_params$knots,
                                density_params$order)
  }

  clr_density <- partial(clr, density_function=density_function)

  if (calculate_norm) {
    norm_clr_density <- integrate(partial(clr_density_squared,
                                          clr_density=clr_density),
                                  lower = 0,
                                  upper = 1)
  }
  else {
    norm_clr_density <- NULL
  }

  return(
    list(density_function = density_function,
         clr_density_function = clr_density,
         norm_clr_density = norm_clr_density
         )
    )
}


# true density to interpolate (we use a beta distribution)
beta_density <- function(quantiles, a, b) {
  dbeta(quantiles, a, b)
}


line_density <- function(quantiles, slope = 0.1) {
  b <- 1 - 0.5 * slope
  return(b + slope * quantiles)
}


truncated_normal_density <- function(quantiles, mean=0, sd=1) {
  dtruncnorm(quantiles, mean, sd)
}


clr <- function(density_function, quantiles) {
  integral_log_density <- integrate(function(x) log(density_function(x)),
                           lower = 0, upper = 1)$value
  log(density_function(quantiles)) - integral_log_density
}


clr_density_squared <- function(clr_density, quantiles) {
  clr_density(quantiles)^2
}


inverse_clr <- function(clr_density, quantiles){
  integral_exp_clr <- integrate(function(x) exp(clr_density(x)),
                                lower = 0, upper = 1)$value

  exp(clr_density(quantiles)) / integral_exp_clr
}

get_spline_clr_density <- function(theta, knots, order) {
  partial(spline_clr_density, theta = theta, knots = knots, order = order)
}


spline_clr_density <- function(quantiles, theta, knots, order) {
  design_matrix <- constrained_spline_design_matrix(x = quantiles, knots = knots, ord = ord)
  design_matrix %*% theta
}


spline_density <- function(quantiles, theta, knots, order) {
  spline_clr_density <- get_spline_clr_density(theta, knots, order)
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


sample_from_density <- function(density_function, n_samples, bins, sample_mode = c("bin", "value")) {
  sample_mode <- match.arg(sample_mode)
  n_bins <- length(bins) - 1
  probs <- get_bin_probabilities(n_bins, density_function)

  if (sample_mode == "value") {
    samples <- rep(NA, n_samples)
    for (i in seq_len(n_samples)) {
      random_bin <- sample(n_bins, 1, prob = probs)
      samples[i] <- runif(1, min = bins[random_bin], max = bins[random_bin + 1])
    }
  }
  else {
    samples <- rmultinom(n = 1, size = N, prob = probs)
  }

  samples

}


