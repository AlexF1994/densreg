smooth.construct.md.smooth.spec <- function (object, data, knots)
{
  # this works for now only for the discrete case where exactly the range is [0,1] and
  # the discrete points are exactly 0 and 1
  # for now I assume the data is sorted by share value
  if (length(object$p.order) == 1)
    m <- rep(object$p.order, 2)
  else m <- object$p.order  # m[1] - basis order, m[2] - penalty order
  m[is.na(m)] <- 2 ## default
  object$p.order <- m
  if (object$bs.dim < 0)
    object$bs.dim <- max(10, m[1] + 1) ## default
  nk <- object$bs.dim - m[1] # basis dimension - order of spline -> number of interior knots
  ## correct dimension of spline to reflect discrete dimensions
  object$bs.dim <- object$bs.dim + 3
  if (nk <= 0)
    stop("basis dimension too small for b-spline order")
  if (length(object$term) != 1)
    stop("Basis only handles 1D smooths")
  x <- data[[object$term]]
  x_cont <- x[x > 0 & x < 1]
  x0 <- x[x == 0]
  x1 <- x[x == 1]
  k <- knots[[object$term]]
  if (is.null(k)) {
    xl <- min(x_cont)
    xu <- max(x_cont)
  }
  else if (length(k) == 2) {
    xl <- min(k)
    xu <- max(k)
    if (xl > min(x) || xu < max(x))
      stop("knot range does not include data")
    if (xl == 0 || xu == 1) {
      stop("Knots 0 and 1 are reserved for discrete
           spline part and shouldn't be provided")
    }
  }
  if (is.null(k) || length(k) == 2) {
    # note that the xl > 0 and xu < 1 at all times and the same must
    # hold for all knots, as 0 and 1 are knots for the discrete part
    # nk = basis dimension - order of spline
    xr <- xu - xl
    xl <- xl / 2
    xu <- xu + (1 - xu) / 2
    kl <- rep(xl, m[1] + 1)
    ku <- rep(xu, m[1] + 1)
    km <- seq(xl ,xu, length = nk)
    k <- c(kl, km, ku)
  }
  else {
    if (length(k) != nk + 2 * m[1] + 2)
      stop(paste("there should be ", nk + 2 * m[1] + 2,
                 " supplied knots"))
    if (0 %in% k || 1  %in% k) {
      stop("Knots 0 and 1 are reserved for discrete
           spline part and shouldn't be provided")
    }
  }
  if (is.null(object$deriv))
    object$deriv <- 0
  # construct design matrix
  design_cont <- splines::spline.des(k, x_cont, m[1] + 2, x_cont * 0 +
                                       object$deriv)$design
  k_discrete <- c(-0.0001, 0 + xl / 2, xu + (1 - xu) / 2 , 1.0001)
  object$knots_discrete <- k_discrete
  design_discrete <- splines::spline.des(k_discrete, x, 1)$design
  n_zeros <- length(x0)
  n_ones <- length(x1)
  matrix_zeros_lower <- matrix(0, ncol = ncol(design_cont),
                               nrow = n_zeros)
  matrix_zeros_upper <- matrix(0, ncol = ncol(design_cont),
                               nrow = n_ones)
  design_cont <- rbind(matrix_zeros_lower, design_cont, matrix_zeros_upper)
  basis_mix <- cbind(design_cont, design_discrete)
  object$X <- basis_mix

  if (!is.null(k)) {
    if (sum(colSums(object$X) == 0) > 0)
      warning("there is *no* information about some basis coefficients")
  }
  if (length(unique(x)) < object$bs.dim)
    warning("basis dimension is larger than number of unique covariates")
  if (is.null(object$mono))
    object$mono <- 0
  if (object$mono != 0) {
    stop("SCOP splines are not supported yet!")
  }
  else {
    cont_dim <- object$bs.dim - 3 # 3 discrete parts
    if (m[2] > 0) {
      S <- diff(diag(cont_dim), differences = m[2])
    }
    else {
      S <- diag(cont_dim)
    }

    matrix_zeros <- matrix(0, ncol = ncol(design_discrete),
                           nrow = cont_dim - m[2])
    S <- cbind(S, matrix_zeros)
    # penalizing the size of the discrete coefficients may be beneficial
    matrix_zeros_bottom <- matrix(0, ncol = ncol(S),
                                  nrow = m[2] + ncol(design_discrete))
    S <- rbind(S, matrix_zeros_bottom)
    ident_matrix <- diag(1, 3)
    S[(nrow(S) - 2): nrow(S), (ncol(S) - 2): ncol(S)] <- ident_matrix
    object$D <- S

    object$S <- list(crossprod(S))
    object$rank <- object$bs.dim - m[2]
    object$null.space.dim <- m[2]
  }
  object$knots <- k
  object$m <- m
  class(object) <- "mdspline.smooth"
  object
}

Predict.matrix.mdspline.smooth <- function (object, data)
{
  m <- object$m[1]
  ll <- object$knots[m + 1]
  ul <- object$knots[length(object$knots) - m]
  m <- m + 1
  x <- data[[object$term]]
  x_cont <- x[x > 0 & x < 1]
  x0 <- x[x == 0]
  x1 <- x[x == 1]
  n_cont <- length(x_cont)
  ind <- x_cont <= ul & x_cont >= ll
  if (is.null(object$deriv))
    object$deriv <- 0
  if (sum(ind) == n_cont) {
    design_cont <- splines::spline.des(object$knots, x_cont, m[1] + 1, x_cont * 0 +
                                         object$deriv)$design
    k_discrete <- object$knots_discrete
    design_discrete <- splines::spline.des(k_discrete, x, 1)$design
    n_zeros <- length(x0)
    n_ones <- length(x1)
    matrix_zeros_lower <- matrix(0, ncol = ncol(design_cont),
                                 nrow = n_zeros)
    matrix_zeros_upper <- matrix(0, ncol = ncol(design_cont),
                                 nrow = n_ones)
    design_cont <- rbind(matrix_zeros_lower, design_cont, matrix_zeros_upper)
    X <- cbind(design_cont, design_discrete)
  }
  else {
    stop("Supplied data is not a density!")
  }
  if (object$mono == 0){
    X
  }
  else {
    stop("SCOP splines are not supported yet!")
  }
}
