# bs = "md" constructs a smoother for the case that the response variable of
# interest is a density in a mixed Bayes Hilbert space. Specification of basis/
# penalization details for the continuous component is as for bs = "ps".
# Specification of the domain and for the discrete component is possible via xt
# (see below).
# Note that the implementation includes the integrate-to-zero constraint of L^2_0,
# i.e., further centering is not reasonable! This means, that the constructor
# should not be used with s() or te(), which always include centering (sum-to-zero)
# constraints for all marginals. Instead, it has to be used with ti, setting
# mc = FALSE for the corresponding component! (It would be desirable to include
# a warning, however, object does not seem to contain information on which kind
# of function it was called in.)
# Special specification via xt-argument, which has to be a list of the same length
# as the vector bs specifying the marginal bases. For all "md" type marginal
# bases, the corresponding xt-list element is again a list with the following
# elements:
# ti(..., xt = list(list(values_discrete = NULL,
#                        weights_discrete = NULL,
#                        domain_continuous = NULL,
#                        penalty_discrete = NULL)))
# - values_discrete: vector of values in the domain of the density that have
#   positive probability mass (dirac measure); defaults to missing (NULL), in
#   which case it is set to c(0, 1); can also be set to be FALSE in which case
#   the discrete component is considered to be empty, i.e., the Lebesgue measure
#   is used as reference measure
# - weights_discrete: vector of weights for the dirac measures corresponding to
#   values_discrete; if missing (NULL) it is set to 1 in all components as default
# - domain_continuous: an interval (i.e., a vector of length 2) specifying the
#   domain of the continuous component of the density; if missing (NULL) it is
#   set to the combined range of values_discrete and the supplied data as default;
#   can also be set to be FALSE in which case the discrete component is considered
#   to be empty, i.e., a sum of dirac measures is used as reference measure
# - penalty_discrete: integer (or NULL) giving the order of differences to be used
#   for the penalty of the discrete component (analogously to m[2] for the continuous
#   component)
# For marginal bases of type "ps"/"bs", the corresponding xt-list element has to
# be set to NULL.

smooth.construct.md.smooth.spec <- function (object, data, knots) {
  if (class(object$xt[[1]]) == "list") {
    if (length(object$xt) == 1) {
      object$xt <- object$xt[[1]]
    } else {
      stop("xt is not specified correctly.")
    }
  }
  discrete <- ifelse(isFALSE(object$xt$values_discrete), FALSE, TRUE)
  continuous <- ifelse(isFALSE(object$xt$domain_continuous), FALSE, TRUE)
  if (!discrete & !continuous) {
    stop("Discrete and continuous components are both empty!")
  }
  x <- data[[object$term]]
  # getting specifications for continuous component
  if (continuous) {
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
  }

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

  # initializing continuous component
  if (!discrete) {
    cont_positions <- seq_along(x)
    x_cont <- x
    if (is.null(object$xt$domain_continuous)) {
      object$xt$domain_continuous <- range(x)
    }
    discrete_dim <- 0
    design_discrete <- Z_discrete <- NULL
  } else if (!continuous) {
    if (any(!(x %in% t_discrete))) {
      stop("Given domain (which is discrete) does not include all observations!")
    }
    discrete_dim <- length(t_discrete)
    x_discrete <- x
    design_cont <- Z_cont <- NULL
    object$bs.dim <- discrete_dim
    cont_dim <- 0
    m <- c(0, 0)
    # cont_positions <- NULL
    # x_cont <- NULL
  } else {
    cont_positions <- which(!(x %in% t_discrete))
    x_cont <- x[cont_positions]
    if (is.null(object$xt$domain_continuous)) {
      object$xt$domain_continuous <- range(t_discrete, x)
    }
  }

  if (continuous) {
    if (object$xt$domain_continuous[1] > min(x_cont) ||
        object$xt$domain_continuous[2] < max(x_cont)) {
      stop("Given domain does not include data corresponding to continuous component.")
    }
  }

  # in mixed case: add discrete value t_{D+1} and weight corresponding to the continuous component
  if (continuous & discrete) {
    t_discrete[length(t_discrete) + 1] <- range(t_discrete, x)[2] + 1
    w_discrete[length(w_discrete) + 1] <- diff(object$xt$domain_continuous)
    discrete_dim <- length(t_discrete)
    # set all continuous values to t_{D+1} for discrete basis
    x_discrete <- x
    x_discrete[cont_positions] <- max(t_discrete)
    object$bs.dim <- object$bs.dim + discrete_dim # combined dimension (before implementing constraints!)
  }

  if (length(unique(x)) < object$bs.dim)
    warning("Basis dimension is larger than number of unique covariates")

  if (continuous) {
    k <- sort(knots[[object$term]])
    if (is.null(k)) {
      xl <- object$xt$domain_continuous[1]
      xu <- object$xt$domain_continuous[2]
    } else if (length(k) == 2) {
      xl <- min(k)
      xu <- max(k)
      if (xl > min(x_cont) || xu < max(x_cont))
        stop("Knot range does not include data corresponding to continuous component.")
    }
    if (is.null(k) || length(k) == 2) {
      # nk = basis dimension - order of spline
      xr <- xu - xl
      # xl <- xl - xr * 0.001 ### E: this shifts the knots slightly outwards, i.e., not (exactly) the supplied knots (for length(k) == 2) are used
      # xu <- xu + xr * 0.001 ### E: I commented it out, since for me it's unclear why not to use the supplied knots
      dx <- (xu - xl)/(nk - 1)
      k <- seq(xl - dx * (m[1] + 1), xu + dx * (m[1] + 1),
               length = nk + 2 * m[1] + 2)
    } else {
      if (length(k) != nk + 2 * m[1] + 2)
        stop(paste("There should be ", nk + 2 * m[1] + 2, " supplied knots"))
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
    }
  }

  if (discrete) {
    k_discrete <- sapply(seq_len(length(t_discrete) + 1),
                         function(j) mean(c(min(t_discrete) - 1, t_discrete, max(t_discrete) + 1)[j:(j+1)]))
    object$knots_discrete <- k_discrete
    object$values_discrete <- t_discrete
    C_discrete <- object$weights_discrete <- w_discrete # integrals of basis functions are equal to weights
    Z_discrete <- MASS::Null(C_discrete)
    design_discrete <- splines::splineDesign(k_discrete, x_discrete, 1)  %*% Z_discrete

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
  }

  object$X <- cbind(design_cont, design_discrete) # combine both design matrices to get final design matrix
  object$Z <- list(cont = Z_cont, discrete = Z_discrete)

  if (!discrete) {
    if (!is.null(k)) {
      if (sum(colSums(object$X) == 0) > 0)
        warning("There is *no* information about some basis coefficients")
    }
  }

  # object$D <- list(continuous = D_cont, discrete = D_discrete)
  # I decided to comment out the definition of object$D since it is not even
  # listed under Value in ?smooth.construct

  if (!discrete) {
    S <- S_cont
    discrete_dim <- 0
    object$xt$penalty_discrete <- NULL
  } else if (!continuous) {
    S <- S_discrete
  } else {
    # combine penalties in a block matrix
    S <- cbind(rbind(S_cont, matrix(0, ncol = ncol(S_cont), nrow = nrow(S_discrete))),
               rbind(matrix(0, nrow = nrow(S_cont), ncol = ncol(S_discrete)), S_discrete))
  }

  object$S <- list(S)
  object$rank <- ifelse(is.null(object$xt$penalty_discrete), cont_dim - m[2],
                        ifelse(object$xt$penalty_discrete == 0,
                               cont_dim - m[2] + discrete_dim - 1, # we loose one dimension when applying constraint
                               cont_dim - m[2] + discrete_dim - object$xt$penalty_discrete))
  object$null.space.dim <- ifelse(is.null(object$xt$penalty_discrete), m[2],
                                  m[2] + object$xt$penalty_discrete)
  if (continuous) {
    object$knots <- k
    object$m <- m
  }
  object$type <- ifelse(!continuous, "discrete",
                        ifelse(!discrete, "continuous", "mixed"))
  class(object) <- "mdspline.smooth"
  return(object)
}

Predict.matrix.mdspline.smooth <- function (object, data) {
  m <- object$m[1] + 2 # object$m[1] + 1
  x <- data[[object$term]]
  t_discrete <- object$values_discrete

  if (object$type != "discrete") {
    ll <- object$xt$domain_continuous[1] # object$knots[m]
    ul <- object$xt$domain_continuous[2] # object$knots[length(object$knots) - m + 1]
  }
  # m <- m + 1
  if (object$type == "continuous") {
    cont_positions <- seq_along(x)
    x_cont <- x
    ind <- list(cont = (x_cont <= ul & x_cont >= ll))
    design_discrete <- NULL
  }
  if (object$type == "mixed") {
    cont_positions <- which(!(x %in% t_discrete[-length(t_discrete)])) # last knot is artificial, corresponding to continuous component
    x_cont <- x[cont_positions]
    x_discrete <- x
    x_discrete[cont_positions] <- max(t_discrete)
    ind <- list(cont = (x_cont <= ul & x_cont >= ll), discrete = (x %in% t_discrete[-length(t_discrete)]))
  }
  if (object$type == "discrete") {
    x_discrete <- x
    ind <- list(discrete = (x %in% t_discrete))
    design_cont <- NULL
  }

  if (is.null(object$deriv)) {
    object$deriv <- 0
  } else if (object$deriv != 0) {
    warning("The mixed density smoother is not intended to be used for derivatives. Reasonable behavior is only guaranteed for deriv = 0.")
  }
  if (sum(sapply(ind, sum)) == length(x)) {
    if (object$type != "discrete") {
      k <- object$knots
      Z_cont <- object$Z$cont
      design_cont <- matrix(0, nrow = length(x), ncol = ncol(Z_cont))
      design_cont[cont_positions, ] <- splines::splineDesign(knots = k, x_cont, ord = m,
                                                             derivs = object$deriv) %*% Z_cont
    }
    if (object$type != "continuous") {
      k_discrete <- object$knots_discrete
      Z_discrete <- object$Z$discrete
      design_discrete <- splines::splineDesign(k_discrete, x_discrete, 1)  %*% Z_discrete
    }
    X <- cbind(design_cont, design_discrete)
  } else {
    stop("Supplied data is not in support of underlying density!")
  }
  if (object$type != "discrete") {
    if (object$mono == 0){
      return(X)
    } else {
      stop("SCOP splines are not supported yet!")
    }
  } else {
    return(X)
  }
}
