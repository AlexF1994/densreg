### The following code contains functions that are used for our analyses

################################################################################
############################# Arithmetic functions #############################
################################################################################

# integrate_num calculates the integral of a vector of function values numerically.
# Arguments:
# x: vector of function values
# type:
# - "joint" assumes, the first and last value of x have positive point mass
#   (dirac measure) and have integration weight 1; the values in between are
#   assumed to be the continuous part (Lebesgue measure), corresponds to the
#   measure mu as in Section 4 of our paper
# - "continuous" integrates with respect to the Lebesgue measure (over all
#   elements in x)
# partial: logical, indicating whether the integral(s) should be calculated over
#          subintervals; If TRUE and type = "joint", the boundary values will be
#          viewed as subintervals
# breaks: if partial = TRUE, this vector specifies the boundary values of the
#         subintervals in (0, 1):
#         (0, breaks[1]); [breaks[1], breaks[2]); ..., (breaks[length(breaks)], 1)
# estimation_points: vector of values, where the functions values (x) are evaluated;
#   Defaults to equidistant values between 0.005 and 0.995
#   The intetgral is calculated using equal weights (diff(estimation_points))

integrate_num <- function(x, type = c("joint", "continuous"), partial = FALSE,
  breaks = NULL, estimation_points = NULL) {
  type <- match.arg(type)
  number_estimation_points <- ifelse(type == "joint", length(x) - 2, length(x))
  if (is.null(estimation_points)) {
    estimation_points <- seq(from = 0.005, to = 0.995, length.out = number_estimation_points)
  }
  stopifnot(length(estimation_points) == number_estimation_points)
  diff_scalar <- diff(estimation_points)[1]
  breaks <- c(0, breaks, 1)
  response <- numeric(length(breaks) + 1)
  if (type == "joint") {
    response[c(1, length(response))] <- x[c(1, length(x))]
    x <- x[- c(1, length(x))]
  }
  for (i in seq_along(breaks[- 1])) {
    interval <- which(estimation_points >= breaks[i] & estimation_points < breaks[i + 1])
    vec <- x[interval]
    diff_vec <- rep(diff_scalar, length(vec))
    response[i + 1] <- sum(vec * diff_vec)
  }
  ifelse(partial,
    ifelse(type == "joint", return(response), return(response[-c(1, length(response))])),
    return(sum(response)))
}

# get_missing_positions computes the column position of the first 7 years (1984-
# 1990), which are missing for East Germany (or regions 5, 6) for each c_age.
# Argument:
#   level: gives the level of spacial refinement ("West_East" for only distinguishing
#          between West and East Germany, "regions", if the finer regions are
#          included)
# (First 33 columns in our matrices correspond to West Germany/region 1, c_age 1,
# next 33 to West Germany/region 1, c_age 2, then West Germany/region 1, c_age 3.
# Then the same for East Germany/region 2 and so on.)
get_missing_positions <- function(level = c("West_East", "regions")) {

  level <- match.arg(level)
  missing_east <- list()

  if (level == "West_East") {
    for(k in 1:3) {
      first_na <- 3 * 33 + (k - 1) * 33 + 1
      last_na <- 3 * 33 + (k - 1) * 33 + 7
      missing_east[[k]] <- first_na:last_na
    }
  } else {
    for (j in 5:6) {
      missing_east[[j]] <- list()
      for(k in 1:3) {
        first_na <- (j - 1) * 3 * 33 + (k - 1) * 33 + 1
        last_na <- (j - 1) * 3 * 33 + (k - 1) * 33 + 7
        missing_east[[j]][[k]] <- first_na:last_na
      }
      missing_east[[j]] <- do.call(c, missing_east[[j]])
    }
  }
  missing_east <- do.call(c, missing_east)
}


# fill_missing_years inserts NA-columns at the specified positions in a matrix
# (in between the existing columns).
# Arguments:
# matrix_without_NA: matrix, where NA columns shall be inserted
# na_positions: vector of positions, where to insert the NA columns
fill_missing_years <- function(matrix_without_NA, na_positions = NULL,
                               level = c("West_East", "regions")) {

  level <- match.arg(level)
  if (is.null(na_positions)) {
    na_positions <- get_missing_positions(level)
  }
  row_number <- nrow(matrix_without_NA)
  col_number <- ncol(matrix_without_NA) + length(na_positions)
  matrix_with_NA <- matrix(rep(NA, row_number * col_number),
                           nrow = row_number, ncol = col_number)
  na_positions <- c(0, na_positions, col_number + 1)
  for (i in seq_along(na_positions[-length(na_positions)])) {
    if (diff(na_positions[c(i, i+1)]) > 1) {
      matrix_with_NA[, (na_positions[i] + 1):(na_positions[i + 1] - 1)] <-
        matrix_without_NA[, (na_positions[i] + 1 - (i - 1)):(na_positions[i + 1] - 1 - (i - 1))]
    } else {
      next
    }
  }
  matrix_with_NA
}

# clr_inv performs the invers clr transformation. (Note: When this code was
# developed, we had not yet implemented the general clr function shown in our
# vignette provided as supplementary material, which can also compute the inverse
# clr transformation).
# Arguments:
# clr_density: matrix containing the clr_densities to transform
# per_col: Are the values of each density stored per column (per_col = TRUE) or
#   per row?
# integrate: Function used to integrate
# ...: Further arguments passed to integrate_num (e.g., estimation_points)
clr_inv <- function(clr_density, per_col = TRUE, integrate = integrate_num, ...) {
  is_vec <- is.vector(clr_density)
  if (is_vec) {
    clr_density <- matrix(clr_density, nrow = ifelse(per_col, length(clr_density), 1))
  }
  density <- exp(clr_density)
  if (per_col) {
    density <- t(density)
  }
  int_density <- apply(density, 1, integrate_num, ...)
  density <- density/int_density
  if (is_vec) {
    return(as.vector(density))
  }
  if (per_col) {
    density <- t(density)
  }
  return(density)
}


# get_discrete_pdf calculates the relative frequencies p_0, p_(0, 1) and p_1
# from a matrix or vector containing mixed densities. Arguments:
# dens_matrix: matrix or vector containing mixed densities
# per_col: Are the values of each density stored per column (per_col = TRUE) or
#   per row? Defaults to FALSE (as is the case for effect matrices)
get_discrete_pdf <- function(dens_matrix, per_col = FALSE) {
  is_vec <- is.vector(dens_matrix)
  if (is_vec) {
    dens_matrix <- matrix(dens_matrix, nrow = ifelse(per_col, length(dens_matrix), 1))
  }
  if (per_col) {
    dens_matrix <- t(dens_matrix)
  }
  dens_matrix_d <- t(apply(dens_matrix, 1, integrate_num, partial = TRUE))
  if (per_col) {
    dens_matrix_d <- t(dens_matrix_d)
  }
  return(dens_matrix_d)
}

################################################################################
################################ Plot functions ################################
################################################################################

# my_minus is a function to add "-"-shaped points to a plot, similar to, e.g.,
# matpoints(..., pch = "-"). The problem with pch = "-" is that it depends
# on the used font, which can be problematic when creating pdfs as the positioning
# gets imprecise (which gets obvious when also adding lines) and the length can
# not be changed (which we need to distinguish between share = 0 and share = 1).
# Arguments:
# x: vector of x-coordinates of the points to plot
# y: vector (of same length as x) or matrix (with ncol same length as x) of
#   y-coordinates of the points to plot
# w: width of the "-"-symbols; Should have length 1 (all rows of y have the same width)
#   or length equal to nrow(y)
# col, lty, lwd: vectors of colors, line types and widths used for the "-"-symbols
# ...: Further arguments passed to segments
my_minus <- function(x, y, w = 0.5, col = NULL, lty = 1, lwd = 1, ...) {
  if (!is.matrix(y)) {
    if (length(x) == 1) {
      y <- matrix(y, ncol = 1)
    } else {
      y <- matrix(y, nrow = 1)
    }
  }
  stopifnot(length(x) == ncol(y))
  stopifnot(length(w) %in% c(1, nrow(y)))
  if (length(w) == 1) {
    w <- rep(w, nrow(y))
  }
  if (length(lty) == 1) {
    lty <- rep(lty, nrow(y))
  }
  if (length(lwd) == 1) {
    lwd <- rep(lwd, nrow(y))
  }
  if (is.null(col)) {
    col <- 1:nrow(y)
  }

  for (i in seq_len(nrow(y))) {
    segments(x0 = x - w[i]/2, x1 = x + w[i]/2, y0 = y[i, ], y1 = y[i, ],
             col = col[i], lty = lty[i], lwd = lwd[i], ...)
  }
}


# plot_densities is used to plot response densities and predictions (Bayes- and
# clr-level) for our model. Arguments:
# dens_matrix: matrix (usually of dimension 102 x 594) with each column containing
#   the density values of one year, region and c_age. The first 33 columns
#   correspond to region 1, c_age 1, next 33 to region 1, c_age 2,
#   then region 1, c_age 3. Then the same for region 2 and so on. The matrix
#   should be named "share" for initial densities and "pred" for predictions,
#   or "..._clr" if it contains clr-transformed densities/predictions
# pdf: boolean value indicating if the plot shall be saved (as pdf)
# name, width, height, path: indicating the name (with the respective path),
#   width and height of the file, if pdf = TRUE;
#   if NULL, name is created automatically out of the name of the matrix given to
#   dens_matrix
# color, ylim, lty, main, ylab: specify the color gradient (each year corresponds
#   to one color), ylim, the line type, the main, and the y-axis limits of the plot
# estimation_points: position of the grid points, where the densities are evaluated
# regions: which regions are supposed to be plotted (defaults to all)
# c_age: which values of c_age are supposed to be plotted (defaults to all)
# cols: which columns are supposed to be plotted (defaults to all)
# single: boolean value indicating whether every year (column) shall be plotted
#   in a single plot or all years in one plot (separately for each region and
#   c_age -> 18 plots in total)
# mod: boolean value; if TRUE "_mod" is added to the file name, otherwise nothing
#   changes in the function (this is to adapt the y-range for examples in paper)
# xaxt: should the x-axis be plotted? (as in plot())
# col_ticks: vector of length 3 specifying the colors of the ticks and labels
#   for s = 0, s in (0, 1), s = 1 on the x-axis; The second color is also used
#   for the axis line between 0 and 1
# mar: as in par()
# ...: Further arguments passed to matplot/plot (e.g. yaxt, las, etc.)
plot_densities <- function(dens_matrix, pdf = TRUE, name = NULL,
  width = 5, height = 3.5, path = "./Images/", color = viridis::plasma(33),
  ylim = NULL, lty = rep(c(1, 2, 4, 5), 9), main = NULL,
  ylab = NULL, estimation_points = seq(from = 0.005, to = 0.995, by = 0.01),
  regions = c("northwest", "west", "southwest", "south", "east", "northeast"),
  c_age = 1:3, cols = 1:33, single = FALSE, mod = FALSE, xaxt = c("s", "n"),
  col_ticks = c("chartreuse2", "mediumblue", "red2"), mar = NULL, cex = 1,
  xlab = "s: income share earned by the woman", ...) {

  xaxt <- match.arg(xaxt)

  # get labels
  matrix_name <- deparse(substitute(dens_matrix))
  clr <- grepl("clr", matrix_name)
  pred <- grepl("pred", matrix_name)
  res <- grepl("res", matrix_name)
  sim <- grepl("sim", matrix_name)
  if (is.null(name)) {
    name <- paste0(ifelse(isTRUE(sim), "sim_", ""),
                   ifelse(isTRUE(clr), "clr_", ""),
                   ifelse(isTRUE(pred), "pred",
                          ifelse(isTRUE(res), "res", "density")))
  }

  if (is.null(ylab)) {
    ylab <- paste0(ifelse(isTRUE(pred), "predicted ", ""),
                   ifelse(isTRUE(clr), "clr-", ""),
                   ifelse(isTRUE(res), "residuals" , "densities"),
                   ifelse(isTRUE(sim), " (Simulation)", ""))
  }

  if(is.null(mar)) {
    mar <- c(5, 4.5, 3, 2) + 0.1
  }

  if (is.null(ylim)) {
    ylim <- c(min(0, min(dens_matrix, na.rm = TRUE)), max(dens_matrix, na.rm = TRUE))
  }

  if (length(col_ticks) == 1) {
    col_ticks <- rep(col_ticks, 3)
  }

  if (any(c("West", "East") %in% regions)) {
    stopifnot("West and East can not be used with specific regions" = length(regions) <= 2)
    all_regions <- c("West", "East")
  } else {
    all_regions <- c("northwest", "west", "southwest", "south", "east", "northeast")
  }
  colnames(dens_matrix) <- rep(NA, ncol(dens_matrix))
  for (j in seq_along(all_regions)) {
    for (k in c_age) {
      colnames(dens_matrix)[((j - 1) * 3 * 33 + (k - 1) * 33 + 1):((j - 1) * 3 * 33 + k * 33)] <-
        paste0(all_regions[j], 1984:2016, "_", k)
    }
  }
  for (j in seq_along(regions)) {
    dens_matrix_mod <- dens_matrix[, which(stringr::str_remove_all(colnames(dens_matrix), "[^(a-zA-Z)]") == regions[j])]
    for (k in c_age) {
      plot_matrix <- dens_matrix_mod[, which(grepl(k, substr(colnames(dens_matrix_mod), nchar(colnames(dens_matrix_mod)), nchar(colnames(dens_matrix_mod)))))]

      # plot
      if (!isTRUE(single)) {
        if (is.null(main)) {
          main_plot <- paste0(#"Region ",
            regions[j], ", cgroup ", k)
        } else {
          main_plot <- main
        }
        if (isTRUE(pdf)) {
          pdf(file = paste0(path, name, "_", regions[j], "_", k, ifelse(isTRUE(mod), "_mod", ""), ".pdf"),
            width = width, height = height)
        }
        par(mar = mar, cex = cex)
        ticks <- seq(0, 1, by = 0.2)
        matplot(estimation_points, plot_matrix[2:101, cols], type = "l", col = color[cols],
          main = main_plot, cex.main = 1, xlab = xlab, ylim = ylim, lty = lty, xlim = c(0, 1),
          ylab = ylab, xaxt = "n", ...)
        if (xaxt == "s") {
          axis(side=1, at = seq(0, 1, by = 0.2), labels = FALSE, lwd.ticks = 0, col = col_ticks[2], ...)
          Map(axis, side=1, at = seq(0.2, 0.8, by = 0.2), col.ticks = col_ticks[2], col.axis = col_ticks[2], ...)
          Map(axis, side=1, at = c(0, 1), labels = format(c(0, 1), nsmall = 1), col.ticks = col_ticks[c(1, 3)], col.axis = col_ticks[c(1, 3)], ...)
        }
        my_minus(c(-0.005, 1.005), t(plot_matrix[c(1, 102), cols]), w = 0.02, col = color[cols], lwd = 0.9)
        abline(h = 0, col = "grey", lwd = 0.5)
        if (isTRUE(pdf)) {
          dev.off()
        }
      } else {
        for (i in cols) {
          if (is.null(main)) {
            main_plot <- paste0("Region ", regions[j], ", cgroup ", k, ", ", i + 1983)
          } else {
            main_plot <- main
          }
          if (isTRUE(pdf)) {
            pdf(file = paste0(path, name, "_", regions[j], "_", k, "_", i + 1983, ifelse(isTRUE(mod), "_mod", ""), ".pdf"),
              width = width, height = height)
          }
          par(mar = mar, cex = cex)
          plot(estimation_points, plot_matrix[2:101, i], type = "l", col = color[i],
            main = main_plot,  cex.main = 1, xlab = xlab,
            ylim = ylim, lty = lty[i], xlim = c(0, 1), ylab = ylab, xaxt = "n", ...)
          if (xaxt == "s") {
            axis(side=1, at = seq(0, 1, by = 0.2), labels = FALSE, lwd.ticks = 0, col = col_ticks[2], ...)
            Map(axis, side=1, at = seq(0.2, 0.8, by = 0.2), col.ticks = col_ticks[2], col.axis = col_ticks[2], ...)
            Map(axis, side=1, at = c(0, 1), labels = format(c(0, 1), nsmall = 1), col.ticks = col_ticks[c(1, 3)], col.axis = col_ticks[c(1, 3)], ...)
          }
          my_minus(c(-0.005, 1.005), t(plot_matrix[c(1, 102), i]), w = 0.02, col = color[i], lwd = 0.9)
          abline( h = 0, col = "grey", lwd = 0.5)
          if (isTRUE(pdf)) {
            dev.off()
          }
        }
      }
    }
  }
}


# plot_histograms_per_year plots relative frequencies for a partition of [0, 1]
# per year on the y-axis (so for one year the values correspond to a histogram).
# Arguments:
# hist_data: Matrix containing the histogram values in the columns (each column corresponds to one year)
# pdf: boolean value indicating if the plot shall be saved (as pdf)
# path: path where to save the file (if pdf = TRUE)
# name: name of the saved file (excluding ".pdf"; if pdf = TRUE)
# width, height: width and height of the pdf (if pdf = TRUE)
# all_years: boolean value indicating if the x-axis should contain all years, i.e.,
#   starting with 1984 (West Germany) or only start with 1991 (East Germany)
# col: color(s) of the different intervals. Should have either length 1 (then the
#   same color is used for all intervals) or number of intervals (nrow(hist_data))
# xcol: color of the ticks and marks of the x-axis; defaults to the color scheme
#   used for the different years in the density plots
# plot_type: indicates whether the data should be plotted as points ("p"), lines ("l")
#   or both ("b"), i.e., points and lines; defaults to points
# breaks: breaks of the histogram intervals in (0, 1); must have length nrow(hist_data) - 2;
#   defaults to NULL corresponding to 0, (0, 1) and 1
# w: If the plot_type includes points and the boundary values should be illustrated
#   as "-" using the function my_minus (see help_functions.R), w is a vector of
#   length 2 giving the width of the "-"-symbol for the s = 0 and s = 1, respectively;
#   If the points for the boundary values should be illustrated as "normal" points,
#   set w = NULL
# mar: as in par()
# pch: point type; If w != NULL this only refers to the points corresponding to
#   the inner intervals; Should have either length 1 (then this point type is used
#   for all (inner) intervals) or number of (inner) intervals
# cex_boundary, cex_inner: point size of the boundary values 0 and 1 and the inner intervals of the histogram intervals
# pch_font: integer giving the font for each histogram "interval"; Should have either length 1 (then the
#   same font is used for all intervals) or number of intervals (nrow(hist_data)); This is important to create
#   "-" of different length to distinguish between the boundary values
# xaxt, yaxt: should the x- and y-axes be plotted? (as in plot())
# y_labels, y_at: Possibility to specify the labels and positions of the y-axis
#   if both are identical (each tick is labeled by the corresponding number), it's
#   sufficient to just give just one of the two arguments
# lty: line type for each histogram "interval"; Should have either length 1 (then the
#   same line type is used for all intervals) or number of intervals (nrow(hist_data))
# legend: boolean value indicating whether a legend should be plotted at the right side of the plot
# ylim: limits of the y-axis; defaults to c(0, max(hist_data))
# text_labels: character containing the text to be added to the plot; Defaults to NULL, i.e., no text
# text_x, text_y: positions of the text, if text_labels != NULL
# ...: Further arguments passed to plot (e.g. main, las, etc.)
plot_histograms_per_year <- function(hist_data, pdf = TRUE, path = NULL, name = NULL,
                                     width = 5, height = 3.5, all_years = TRUE, col = 1,
                                     xcol = viridis::plasma(33), plot_type = c("p", "l", "b"),
                                     breaks = NULL, w = c(0.5, 0.3), mar = NULL,
                                     pch = 1, cex_boundary = 1, cex_inner = 0.7,
                                     xaxt = c("s", "n"), yaxt = c("s", "n"), y_labels = NULL,
                                     y_at = NULL, lty = 1, legend = FALSE, ylim = NULL,
                                     text_labels = NULL, text_x = 0.5, text_y = 0.5,
                                     cex = 1, xlab = "year", ...) {

  plot_type <- match.arg(plot_type)
  year <- seq(ifelse(all_years, 1984, 1991), 2016)
  year_5 <- seq(ifelse(all_years, 1985, 1995), 2015, by = 5)
  xcol_5 <- xcol[seq(ifelse(all_years, 2, 12), 32, by = 5)] # yields colors corresponding to 1985, 1990, ..., 2015 or 1995, 2000, ..., 2015
  n <- nrow(hist_data)

  if (length(pch) == 1) {
    pch_length <- ifelse(is.null(w), n, n - 2)
    pch <- rep(pch, pch_length)
  }
  if (length(lty) == 1) {
    lty <- rep(lty, n)
  }
  if (length(col) == 1) {
    col <- rep(col, n)
  }
  if (is.null(ylim)) {
    ylim <- c(0, max(hist_data, na.rm = TRUE))
  }
  if (!is.null(y_labels) | !is.null(y_at)) {
    yaxt_intern <- "n"
  } else {
    yaxt_intern <- yaxt
  }

  if (pdf) {
    pdf(file = paste0(path, name, ".pdf"), width = width, height = height)
  }
  if (isTRUE(legend) & is.null(mar)) {
    mar <- c(5, 4.5, 3, 10) + 0.1
  } else if (isFALSE(legend) & is.null(mar)) {
    mar <- c(5, 4.5, 3, 2) + 0.1
  }
  par(mar = mar, cex = cex)

  type <- ifelse(plot_type %in% c("p", "l"), plot_type, "p")
  vec <- seq_len(n)[-2]
  plot(year, hist_data[2, ], type = type, col = col[2], xaxt = "n", xlab = xlab,
       ylim = ylim, pch = ifelse(is.null(w), pch[2], pch[1]), cex = cex_inner,
       cex.main = 1, lty = lty[2], yaxt = yaxt_intern, ...)

  if (type == "p") {
    if (!is.null(w)) {
      my_minus(x = year, y = hist_data[c(1, n), ], w = w, col = col[c(1, n)])
      for (i in vec[-c(1, length(vec))]) {
        points(year, hist_data[i, ], pch = pch[i-1], col = col[i], cex = cex_inner)
      }
    } else {
      for (i in vec) {
        points(year, hist_data[i, ], pch = pch[i], col = col[i], cex = cex_inner)
      }
    }
    if (plot_type == "b") {
      for (i in 1:n) {
        lines(year, hist_data[i, ], col = col[i], lty = lty[i], ...)
      }
    }
  } else {
    for (i in vec) {
      lines(year, hist_data[i, ], col = col[i], lty = lty[i], ...)
    }
  }
  abline( h = 0, col = "grey", lwd = 0.5)

  if (xaxt == "s") {
    Map(axis, side=1, at = year, labels = FALSE, tick = TRUE,
        col.ticks = xcol[seq(ifelse(all_years, 1, 8), 33, by = 1)])
    Map(axis, side=1, at = year_5, labels = year_5, lwd.ticks = 2,
      col.ticks = xcol_5, col.axis =  xcol_5)
  }
  if (!is.null(y_labels) & is.null(y_at)) {
    y_at <- y_labels
  } else if(is.null(y_labels) & !is.null(y_at)) {
    y_labels <- y_at
  }
  if (!is.null(y_labels) & !is.null(y_at)) {
    axis(2, at = y_at, labels = y_labels, ...)
  }
  if (isTRUE(legend)) {
    breaks <- c(0, breaks, 1)
    legend <- character(length(breaks) + 1)
    legend[c(1, length(legend))] <- c("0", "1")
    for (i in seq_along(breaks[- 1])) {
      legend[i + 1] <- ifelse(breaks[i] == 0,
                              paste0("(0, ", round(breaks[i + 1], digits = 3), ")"),
                              paste0("[", round(breaks[i], digits = 3), ", ",
                                     round(breaks[i + 1], digits = 3), ")"))
    }
    legend(2018, 0.32, xpd = TRUE, legend = legend, title = "woman's income share s:",
           title.col = 1, text.col = col, lty = 1, pch = pch, bty = "n", col = col)
  }
  if(!is.null(text_labels)) {
    text(x = text_x, y = text_y, labels = text_labels)
  }
  if (pdf) {
    dev.off()
  }
}

# plot_for_region_c_age calls the plot_histograms_per_year function for the
# specific categories, i.e., each combination of region and c_age.
# Arguments:
# j, k: Integers giving the number of the region and c_age, respectively
# all_years, breaks, pdf, path, col, w, pch, cex_boundary, cex_inner, ylim: as in plot_histograms_per_year
# name_add: suffix added to the name of the file
# ...: Further arguments passed to plot_histograms_per_year
plot_for_region_c_age <- function(j, k, all_years = TRUE, breaks = NULL, pdf = TRUE, path = NULL,
                                        col, w = c(0.5, 0.3), pch = 1, cex_boundary = 1,
                                        cex_inner = 0.7, ylim = NULL, name_add = "_hist_years", ...) {
  ylab_initial <- "relative frequency"
  ylab_pred <- "predicted relative frequency"
  intercept_West_East <- ifelse(all_years, 1, ifelse(j < 5,  1, 8))
  first_year <- (j - 1) * 3 * 33 + (k - 1) * 33 + intercept_West_East
  last_year <- (j - 1) * 3 * 33 + k * 33
  # Initial densities
  initial_k <- apply(share[, first_year:last_year], 2, integrate_num, partial = TRUE,
                       breaks = breaks)
  # Predictions
  pred_k <- apply(pred[, first_year:last_year], 2, integrate_num, partial = TRUE, breaks = breaks)

  ### Plot histograms
  if (is.null(ylim)) {
    ylim <- c(0, max(c(initial_k, pred_k), na.rm= TRUE))
  }

  # regions <- c("north-west", "nrw", "west", "south", "east", "north-east"),
  regions <- c("northwest", "west", "southwest", "south", "east", "northeast")
  main_k <- paste0("Region ", regions[j], ", cgroup ", k)
  # Initial densities
  plot_histograms_per_year(hist_data = initial_k, path = path, breaks = breaks,
                           name = paste0("initial_", regions[j], "_", k, name_add),
                           col = col, main = main_k, ylab = ylab_initial,
                           pch = pch, plot_type = "p", all_years = all_years, pdf = pdf,
                           cex_boundary = cex_boundary, cex_inner = cex_inner, ylim = ylim, ...)
  # Predictions
  plot_histograms_per_year(hist_data = pred_k, path = path, breaks = breaks,
                           name = paste0("pred_", regions[j], "_", k, name_add),
                           col = col, main = main_k, ylab = ylab_pred,
                           all_years = all_years, plot_type = "l", pdf = pdf, ylim = ylim, ...)
  invisible(0)
}

# plot_effects plots an estimated model effect. Arguments:
# effect_matrix: matrix containing the estimated effect
# pdf: boolean value indicating if the plot shall be saved (as pdf)
# name, width, height, path: indicating the name, width and height of the file,
#   if pdf = TRUE;
#   if name = NULL, the name is created automatically out of the name of the
#   matrix given to dens_matrix
# type: either "joint" or "continuous" specifying whether the densities are with
#   respect to the mixed reference measure or the Lebesgue measure
# effect_matrix_d: matrix containing the effects of the discrete model; If passed,
#   all discrete values (including the one at 0.5 corresponding to dual-earners)
#   are added to the plot as "-". Defaults to NULL, i.e., only boundary is plotted
# pch, my_minus_lwd, my_minus_lty, cex_inner: Arguments for point specification of
#   dual-earner discrete value; if pch = "minus", the my_minus function is used
#   and cex_inner is ignored (default); my_minus_lwd and my_minus_lty specify the
#   width and the linetype of the minus symbol
# legend: containing the legend argument passed to function legend. If missing (default),
#   no legend is added
# legend_x, legend_y: x and y-co-ordinates to be used to position the legend.
# legend_y.intersp: character interspacing factor for vertical (y) spacing in legend.
# legend_ncol: the number of columns in which to set the legend items.
# legend_lty: line types for lines appearing in the legends. Defaults to NULL, i.e.,
#   no lines in legend.
# legend_title: title to be placed at the top of the legend. Defaults to no title.
# legend_cex: cex to be used in legend
# legend_col: defaults to col
# legend_seg.len: the length of lines drawn to illustrate lty and/or lwd in legend
# main: contains the text for the plot title; If null, no title is added
# col, lty, ylab, ylim: specify the color vector, the line type and the label and
#   range on the y axis of the plot
# mar: graphic parameter as in par(). If NULL (default), mar is chosen depending
#   on the existence of a main and the value of odds
# odds: If true, the margin to the left is extended by two to leave enough space
#   for the fraction in the label of the y-axis
# G: number of bins for continuous component
# abline_h, abline_v, abline_col, abline_lty: values, where to draw horizontal/vertical
#   lines and their colors and line types
# add_geom_mean, lty_geom_mean: should the geometric mean of the densities be calculated
#   and added to the plot as horizontal lines? (Only reasonable for Bayes level)
# yaxt: should the y axis been drawn? ("s": yes, "n": no)
# y_labels, y_at: Possibility to specify the labels and positions of the y-axis
#   if both are identical (each tick is labeled by the corresponding number), it's
#   sufficient to just give just one of the two arguments
# ...: Further plot arguments (passed to matplot)
plot_effects <- function(effect_matrix, pdf = TRUE, name = NULL, width = 5,
  height = 3, path = "./Images/", type = c("joint", "continuous"), effect_matrix_d = NULL,
  pch = "minus", radius = 0.01, my_minus_lwd = NULL, my_minus_lty = NULL, my_minus_w = 0.02,
  cex_inner = 0.7, legend = NULL,
  legend_x = NULL, legend_y = NULL, legend_y.intersp = 1, legend_x.intersp = 1,
  legend_ncol = 1, legend_lty = NULL, legend_title = NULL, legend_cex = 1, legend_col = NULL,
  legend_lwd = integer(0),
  legend_seg.len = 2, main = NULL, col, lty = 1, lwd = 1, ylab = NULL, ylim = NULL,
  xlab = "s: income share earned by the woman", mar = NULL,
  odds = FALSE, G = 100, abline_h = 0, abline_v = NULL, abline_col = "grey", abline_lty = 1,
  add_geom_mean = FALSE, lty_geom_mean = 1, yaxt = "s", y_labels = NULL, y_at = NULL, ...) {

  if(!is.matrix(effect_matrix)) {
    effect_matrix <- matrix(effect_matrix, nrow = 1)
  }

  type <- match.arg(type)

  step_size <- 1 / G
  estimation_points = seq(step_size / 2, 1 - step_size / 2, length.out = G)

  # get labels
  matrix_name <- deparse(substitute(effect_matrix))
  clr <- !(grepl("pdf", matrix_name))
  sim <- grepl("sim", matrix_name)
  if (is.null(name)) {
    name <- paste0("estimated_", sub("_pdf", "", sub(".*\\$", "", matrix_name)),
                   ifelse(isTRUE(clr), "_clr", ""), ifelse(isTRUE(sim), "_sim", ""))
  }

  if (is.null(mar)) {
    if (is.null(main)) {
      mar <- c(5, 5.3, 0.5, 1) + 0.1
    } else {
      mar <- c(5, 5.3, 4, 3) + 0.1
    }
    if (odds) {
      mar <- mar + c(0, 2, 0, 0)
    }
  }

  if (length(lty) != nrow(effect_matrix)) {
    lty <- rep(lty, nrow(effect_matrix))
  }
  if (length(lwd) != nrow(effect_matrix)) {
    lwd <- rep(lwd, nrow(effect_matrix))
  }

  if (is.null(ylab)) {
    ylab <- paste0("estimated ", stringr::str_replace_all(name, "_", " "))
  }

  if (is.null(ylim)) {
    if (isTRUE(clr)) {
      ylim <- range(effect_matrix, na.rm = TRUE)
    } else {
      ylim <- c(0, max(effect_matrix, na.rm = TRUE))
    }
  }

  if (!is.null(y_labels) | !is.null(y_at)) {
    yaxt_intern <- "n"
  } else {
    yaxt_intern <- yaxt
  }

  if (type == "joint") {
    continuous <- 2:(ncol(effect_matrix) - 1)
  } else {
    continuous <- 1:ncol(effect_matrix)
  }

  if (nrow(effect_matrix) == 1) {
    plot_matrix <- matrix(effect_matrix[, continuous], nrow = length(continuous))
  } else {
    plot_matrix <- t(effect_matrix[, continuous])
  }

  if (isTRUE(pdf)) {
    pdf(paste0(path, name, ".pdf"), width = width, height = height)
  }
  par(mar = mar)
  matplot(estimation_points, plot_matrix, type = "l", main = main, lty = lty, lwd = lwd,
          ylim = ylim, ylab = ylab, xlab = xlab, col = col, yaxt = yaxt_intern, ...)
  abline(h = abline_h, v = abline_v, col = abline_col, lty = abline_lty, lwd = 0.5)

  if (!is.null(y_labels) & is.null(y_at)) {
    y_at <- y_labels
  } else if(is.null(y_labels) & !is.null(y_at)) {
    y_labels <- y_at
  }
  if (!is.null(y_labels) & !is.null(y_at)) {
    axis(2, at = y_at, labels = y_labels, ...)
  }

  if (add_geom_mean) {
    if (clr) {
      warning("Geometric mean not reasonable on clr-level!")
    }
    if (type == "joint") {
      geom_mean <- apply(effect_matrix, 1, function(x) exp(1/3 * integrate_num(log(x))))
    } else {
      geom_mean <- apply(effect_matrix, 1, function(x) exp(integrate_num(log(x), type = "continuous")))
    }

    abline(h = geom_mean, col = col, lty = lty_geom_mean, lwd = 0.5)
  }
  if (is.null(my_minus_lty)) {
    my_minus_lty <- 1
  }
  if (is.null(my_minus_lwd)) {
    my_minus_lwd <- lwd
  }
  if (type == "joint") {

    my_minus(c(-0.005, 1.005), effect_matrix[, c(1, ncol(effect_matrix))], w = my_minus_w,
             col = col, lwd = my_minus_lwd, lty = my_minus_lty)
  }
  if (!is.null(effect_matrix_d)) {
    if (pch == "minus") {
      my_minus(0.5, effect_matrix_d[, 2], w = my_minus_w, col = col, lwd = my_minus_lwd, lty = my_minus_lty)
    } else if (pch == "circle") {
      sapply(seq_len(nrow(effect_matrix_d)),
             FUN = function(i) plotrix::draw.circle(0.5, effect_matrix_d[i, 2], radius = radius,
                                                    border = col[i], lty = lty[i], lwd = lwd[i]))
    } else {
      points(rep(0.5, nrow(effect_matrix_d)), effect_matrix_d[, 2], col = col, pch = pch, cex = cex_inner)
    }
  }

  if (!is.null(legend)) {
    if (is.null(legend_col)) {
      legend_col <- col
    }
    if (!(is.null(legend_lty)) & length(legend_lwd) == 0) {
      legend_lwd <- 1
    }
    if (is.null(legend_lwd)) {
      legend_lwd <- lwd
    }
    legend(x = legend_x, y = legend_y, legend = legend, text.col = legend_col,
           cex = legend_cex, bty = "n",
           y.intersp = legend_y.intersp,
           x.intersp = legend_x.intersp, ncol = legend_ncol,
           xpd = TRUE, seg.len = legend_seg.len,
           title = legend_title, title.col = "black", title.adj = 0)
  }

  if (isTRUE(pdf)) {
    dev.off()
  }

  if (add_geom_mean) {
    return(geom_mean)
  }
}

# multiply_borders is a helper for the mixed heatmaps that copies the borders,
# i.e., the first and last rows and columns of the matrix to plot to make them
# appear wider in the resulting plot.
# Arguments:
# odds_matrix: matrix containing the values for the heatmap to plot
# n_border_left_right: value how often to copy the first/last column
# n_border_up_down: value how often to copy the first/last row
multiply_borders <- function(odds_matrix, n_border_left_right, n_border_up_down) {
  n_col <- ncol(odds_matrix)
  n_row <- nrow(odds_matrix)
  odds_matrix_multiplied <- cbind(matrix(odds_matrix[, 1], nrow = n_row,
                                         ncol = n_border_left_right),
                                  odds_matrix[, 2:(n_col - 1)],
                                  matrix(odds_matrix[, n_col], nrow = n_row,
                                         ncol = n_border_left_right))
  odds_matrix_multiplied <- rbind(matrix(odds_matrix_multiplied[1, ], nrow = n_border_up_down,
                                         ncol = ncol(odds_matrix_multiplied), byrow = TRUE),
                                  odds_matrix_multiplied[2:(n_col - 1), ],
                                  matrix(odds_matrix_multiplied[n_col, ], nrow = n_border_up_down,
                                         ncol = ncol(odds_matrix_multiplied), byrow = TRUE))
  return(odds_matrix_multiplied)
}

# get_color is a helper for the mixed heatmaps that returns the corresponding color
# of a value in a heatmap. Arguments:
# x: values (matrix or vector) to determine color for
# image_col: vector of colors to use
# image_breaks: vector of values, where the colors change (including boundary
#   values); must have length(image_col) + 1
get_color <- function(x, image_col, image_breaks) {
  get_pos <- function(k, x) { image_breaks[k] <= x & image_breaks[k + 1] > x }
  if (is.vector(x) & length(x) == 1) {
    pos <- which(sapply(seq_along(image_breaks), get_pos, x = x))
  } else {
    pos <- apply(sapply(seq_along(image_breaks), get_pos, x = x), 1, which)
  }
  if (is.matrix(x)) {
    return(matrix(image_col[pos], nrow = nrow(x), ncol = ncol(x)))
  } else {
    return(image_col[pos])
  }
}

# plot_odds plots the mixed heatmaps used for the illustration of the (log) odds.
# Arguments:
# odds_matrix: matrix containing the values of the heatmap to plot
# odds_matrix_d: discrete component heatmap (optional)
# mar, main, xlab, ylab: as in function plot()
# x_axis, y_axis: boolean indicating whether to draw an x-/y-axis
# n_border_left_right, n_border_up_down: as in function multiply_borders()
# share: vector containing the positions of the grid points (x- and y-axis), where
#   the heatmap is evaluated
# lwd_seg: width of dotted lines to separate outer from inner bands
# image_col, image_breaks: as col and breaks in function image()
# contour_cex, image_breaks_contour: values for cex and levels to pass to contour
# axis_cex: cex value to pass to axis()
# legend: containing the legend argument passed to function legend() (usually used
#   for c_age). If missing (default), no legend is added
# legend_x, legend_y: x and y-co-ordinates to be used to position the legend.
# legend_year: containing the legend argument passed to a second function legend()
#   (usually used for year). If missing (default), no legend is added
# legend_year_x, legend_year_y: x and y-co-ordinates to be used to position the
#   second legend.
# legend_cex: cex value used for both legends

plot_odds <- function(odds_matrix, odds_matrix_d = NULL, mar, main = "", xlab = "",
                      ylab = "", x_axis = FALSE, y_axis = FALSE, n_border_left_right,
                      n_border_up_down, share = seq(from = 0.005, to = 0.995, by = 0.01),
                      lwd_seg = 0.8, image_col, image_breaks, contour_cex = 1,
                      image_breaks_contour = image_breaks, axis_cex = 1,
                      legend = NULL, legend_x = 1.1, legend_y = 0.5, legend_year = NULL,
                      legend_year_x = 0.5, legend_year_y = 1.15, legend_cex = 1) {
  share_diff <- diff(share)[1]
  left <- share[1] - rev(seq_len(n_border_left_right)) * share_diff
  n <- length(share)
  right <- share[n] + seq_len(n_border_left_right) * share_diff
  bottom <- share[1] - rev(seq_len(n_border_up_down)) * share_diff
  top <- share[n] + seq_len(n_border_up_down) * share_diff

  if (!is.null(odds_matrix_d)) {
    width_rect <- n_border_left_right * share_diff
    height_rect <- n_border_up_down * share_diff
    xlim <- c(-2 * width_rect, 1 + 2 * width_rect)
    ylim <- c(-2 * height_rect, 1 + 2 * height_rect)
  } else {
    xlim <- c(left[1], right[n_border_left_right])
    ylim <- c(bottom[1], top[n_border_up_down])
  }

  odds_matrix_mod <- multiply_borders(odds_matrix, n_border_left_right, n_border_up_down)
  par (mar = mar)
  image(c(left, share, right), c(bottom, share, top), t(odds_matrix_mod),
        col = image_col, main = main, xaxt = "n", yaxt = "n", xlim = xlim,
        ylim = ylim, xlab = xlab, ylab = ylab, breaks = image_breaks)
  contour(share, share, t(odds_matrix[2:101, 2:101]), levels = image_breaks_contour,
          add = TRUE, labcex = contour_cex)
  if (!is.null(odds_matrix_d)) {
    col <- get_color(odds_matrix_d, image_breaks = image_breaks,
                     image_col = image_col)
    col <- col[c(1:3, 7:9 , 1, 4, 7, 3, 6, 9)]
    rect(xleft = c(rep(-2 * width_rect, 3), # left ribbon
                   rep(1 + width_rect, 3), # right ribbon
                   rep(c(-2 * width_rect, 0, 1), 2)), # bottom top ribbon and bottom ribbon
         xright = c(rep(-width_rect, 3), # left ribbon
                    rep(1 + 2 * width_rect, 3), # right ribbon
                    rep(c(0, 1, 1 + 2 * width_rect), 2)), # bottom top ribbon and bottom ribbon
         ybottom = c(rep(c(- height_rect, 0, 1), 2), # left and right ribbon
                     rep(-2 * height_rect, 3), # bottom ribbon
                     rep(1 + height_rect, 3)), # top ribbon
         ytop = c(rep(c(0, 1, 1 + height_rect), 2),# left and right ribbon
                  rep(- height_rect, 3), # bottom ribbon
                  rep(1 + 2 * height_rect, 3)), # top ribbon
         col = col, border = NA)
    segments(x0 = c(-width_rect, 1 + width_rect, rep(0, 2)),
             y0 = c(rep(0, 2), -height_rect, 1 + height_rect),
             x1 = c(-width_rect, 1 + width_rect, rep(1, 2)),
             y1 = c(rep(1, 2), -height_rect, 1 + height_rect),
             lty = 3, lwd = lwd_seg)
  }
  abline(h = c(0, 1), v = c(0, 1))
  if (x_axis) {
    axis(1, at = seq(0, 1, by = 0.1), lwd = 0, lwd.ticks = 1, cex.axis = axis_cex)
  }
  if (y_axis) {
    axis(2, at = seq(0, 1, by = 0.1), lwd = 0, lwd.ticks = 1, las = 2, cex.axis = axis_cex)
  }
  if (!is.null(legend)) {
    legend(x = legend_x, y = legend_y, legend = legend, xpd = TRUE, bty = "n",
           xjust = 0.5, yjust = 0.5, cex = legend_cex)
  }
  if (!is.null(legend_year)) {
    legend(x = legend_year_x, y = legend_year_y, legend = legend_year,
           xpd = TRUE, bty = "n", xjust = 0.5, yjust = 0.5, cex = legend_cex)
  }
}

################################################################################
######################## Functions to get model effects ########################
################################################################################

# get_single_effect computes the model effect matrix from a specific effect from
# an the predicted terms of a model estimated with mgcv::gam
# Arguments:
# pred_terms: object created via predict(..., type = "terms")
# positions: column positions of the effects of interest in pred_terms
# G: number of discrete values + number of bins for continuous component
get_single_effect <- function(pred_terms, positions, G = 102) {
  t(unique(Reduce("+", lapply(positions,
                            function(k) matrix(pred_terms[, k], nrow = G))),
         MARGIN = 2))
}

# get_estimated_model_effects calculates the model effects (Bayes- and clr-level).
# Arguments:
# model: gam model object
# G: number of bins for continuous component
get_estimated_model_effects <- function(model, G = 100) {
  # Package to fit functional regression models
  library(mgcv)
  # # Package to construct a matrix of different class (we need it to calculate the
  # # predictions for the region effect)
  # library(Matrix)

  ###################### Get estimated effects from model ######################
  # predict(..., type = "terms") contains estimated effects on clr-level with columns:
  # - 1: intercepts per covariate combination (not of interest for us)
  # - 2: smooth intercept (reference West Germany, c_age other)
  # - 3: smooth effect for East Germany
  # - 4-5: smooth effects for c_age 2 (7-18) and 1 (0-6)
  # - 6-10: smooth interaction of West/East Germany and c_age (3)/2/1
  #         (3 only for East Germany!)
  # - 11: smooth effect for year
  # - 12: smooth interaction of year and East Germany
  # - 13-14: smooth interaction of year and c_age 2/1
  # - 15-19: smooth interaction of year and West/East Germany and c_age (3)/2/1
  #          (3 only for East Germany!)
  pred_terms <- predict(model, type = "terms")
  # share <- seq(step_size / 2, 1 - step_size / 2, by = step_size)
  # G <- length(share) + 2
  G <- G + 2
  stopifnot("Given number of discrete values + number of bins for continuous component is not compatible with given model" = nrow(pred_terms) %% G == 0)

  ### Intercept (beta_0)
  intercept_origcod <- c(get_single_effect(pred_terms, 2, G))

  ### West/East (beta_{West_East})
  West_East_origcod <- get_single_effect(pred_terms, 3, G)

  ### Region (beta_{region})
  # region_origcod <-

  ### c_age (beta_{c_age})
  c_age_origcod <- get_single_effect(pred_terms, 4:5, G)[3:1,] # change order of c_age to 1, 2, 3

  ### Interaction c_age/West_East (beta_{c_age, West_East})
  West_East_c_age_origcod <- get_single_effect(pred_terms, 6:10, G)[c(3:1, 6:4),] # change order of c_age to 1, 2, 3 for West and East

  ### Year intercept (g_0(year))
  year_intercept_origcod <- get_single_effect(pred_terms, 11, G)

  ### Year West_East (g_{West_East}(year))
  year_West_origcod <- matrix(0, ncol = G, nrow = 33) # reference
  year_East_origcod <- get_single_effect(pred_terms, 12, G)[-1,] # first row contains only zeros (for West Germany)
  # add NA for years 1984-1990 (for compability with number of columns in West)
  year_East_origcod <- rbind(matrix(NA, ncol = G, nrow = 7), year_East_origcod)

  ### Year c_age (g_{c_age}(year))
  # Year, c_age 3 (g_{3}(year))
  year_3_origcod <- matrix(0, ncol = G, nrow = 33)  # reference
  # Year, c_age 2 (g_{2}(year))
  year_2_origcod <- get_single_effect(pred_terms, 13, G)[-1,] # first row contains only zeros (for West Germany)
  # Year, c_age 1 (g_{1}(year))
  year_1_origcod <- get_single_effect(pred_terms, 14, G)[-1,] # first row contains only zeros (for West Germany)

  ### Interaction c_age/West_East/year (g_{c_age, West_East}(year))
  # Year, West, c_age 3 (g_{3, West}(year))
  year_West_3_origcod <- matrix(0, ncol = G, nrow = 33)  # reference
  # Year, West, c_age 2 (g_{2, West}(year))
  year_West_2_origcod <- get_single_effect(pred_terms, 15, G)[-1,] # first row contains only zeros (for West Germany)
  # Year, West, c_age 1 (g_{1, West}(year))
  year_West_1_origcod <- get_single_effect(pred_terms, 16, G)[-1,] # first row contains only zeros (for West Germany)
  # Year, East, c_age 3 (g_{3, East}(year))
  year_East_3_origcod <- get_single_effect(pred_terms, 17, G)[-1,] # first row contains only zeros (for West Germany)
  # add NA for years 1984-1990 (for compability with number of columns in West)
  year_East_3_origcod <- rbind(matrix(NA, ncol = G, nrow = 7), year_East_3_origcod)
  # Year, East, c_age 2 (g_{2, East}(year))
  year_East_2_origcod <- get_single_effect(pred_terms, 18, G)[-1,] # first row contains only zeros (for West Germany)
  # add NA for years 1984-1990 (for compability with number of columns in West)
  year_East_2_origcod <- rbind(matrix(NA, ncol = G, nrow = 7), year_East_2_origcod)
  # Year, East, c_age 1 (g_{1, East}(year))
  year_East_1_origcod <- get_single_effect(pred_terms, 19, G)[-1,] # first row contains only zeros (for West Germany)
  # add NA for years 1984-1990 (for compability with number of columns in West)
  year_East_1_origcod <- rbind(matrix(NA, ncol = G, nrow = 7), year_East_1_origcod)


  ############# Transform to reference coding used in first paper ##############
  # Reference category: West, 3, 1991

  ### intercept = beta_0 + g_0(1991)
  intercept <- intercept_origcod + year_intercept_origcod[8, ]

  ### West_East = beta_{West_East} + beta{3, West_East}
  ###           + g_{West_East}(1991) + g_{3, West_East}(1991)
  West_East <- West_East_origcod + West_East_c_age_origcod[c(3, 6), ] +
    matrix(c(year_West_origcod[8, ], year_East_origcod[8, ]), byrow = TRUE, nrow = 2) +
    matrix(c(year_West_3_origcod[8, ], year_East_3_origcod[8, ]), byrow = TRUE, nrow = 2)

  # ### region = beta_{region}
  # region <- region_origcod

  ### c_age = beta_{c_age} + beta_{c_age, West}
  ###         + g_{c_age}(1991) + g_{c_age, West}(1991)
  c_age <- c_age_origcod + West_East_c_age_origcod[1:3, ] +
    matrix(c(year_1_origcod[8, ], year_2_origcod[8, ], year_3_origcod[8, ]),
           byrow = TRUE, nrow = 3) +
    matrix(c(year_West_1_origcod[8, ], year_West_2_origcod[8, ], year_West_3_origcod[8, ]),
           byrow = TRUE, nrow = 3)

  ### c_age_West_East = beta_{c_age, West_East} - beta_{3, West_East} - beta_{c_age, West}
  ###                 + g_{c_age, West_East}(1991) - g_{3, West_East}(1991) - g_{c_age, West}(1991)
  West_East_c_age <- West_East_c_age_origcod -
    matrix(c(rep(West_East_c_age_origcod[3, ], 3),
             rep(West_East_c_age_origcod[6, ], 3)), byrow = TRUE, nrow = 6) -
    matrix(rep(t(West_East_c_age_origcod[1:3, ]), 2), byrow = TRUE, nrow = 6) +
    matrix(c(year_West_1_origcod[8, ], year_West_2_origcod[8, ], year_West_3_origcod[8, ],
             year_East_1_origcod[8, ], year_East_2_origcod[8, ], year_East_3_origcod[8, ]),
           byrow = TRUE, nrow = 6) -
    matrix(c(rep(year_West_3_origcod[8, ], 3), rep(year_East_3_origcod[8, ], 3)),
           byrow = TRUE, nrow = 6) -
    matrix(rep(c(year_West_1_origcod[8, ], year_West_2_origcod[8, ], year_West_3_origcod[8, ]), 2),
           byrow = TRUE, nrow = 6)

  ### year = g_0(year) - g_0(1991)
  year_intercept <- t(apply(year_intercept_origcod, 1, "-",
                            year_intercept_origcod[8, ]))

  ### year_{West_East} = g_{West_East}(year) + g_{3, West_East}(year)
  ###                  - g_{West_East}(1991) - g_{3, West_East}(1991)
  # g_{West}(year) is constant 0 for all years; already in original coding
  year_West <- year_West_origcod
  # g_{East}(year)
  year_East <- t(apply(year_East_origcod, 1, "-",
                       year_East_origcod[8, ] + year_East_3_origcod[8, ])) +
    year_East_3_origcod

  ### year_{c_age} = g_{c_age}(year) + g_{c_age, West}(year)
  ###                - g_{c_age}(1991) - g_{c_age, West}(1991)
  # g_1(year)
  year_1 <- t(apply(year_1_origcod, 1, "-",
                    year_1_origcod[8, ] + year_West_1_origcod[8, ])) + year_West_1_origcod
  # g_2(year)
  year_2 <- t(apply(year_2_origcod, 1, "-",
                    year_2_origcod[8, ] + year_West_2_origcod[8, ])) + year_West_2_origcod
  # g_3(year) is constant 0 for all years; already in original coding
  year_3 <- year_3_origcod

  ### year_{c_age, West_East} = g_{c_age, West_East}(year) - g_{3, West_East}(year)
  ###                         - g_{c_age, West}(year)
  ###                         - g_{c_age, West_East}(1991) + g_{3, West_East}(1991)
  ###                         + g_{c_age, West}(1991)
  # g_{1, West}(year) is constant 0 for all years
  year_West_1 <- matrix(0, nrow = nrow(year_West_1_origcod), ncol = ncol(year_West_1_origcod))
  # g_{2, West}(year) is constant 0 for all years
  year_West_2 <- matrix(0, nrow = nrow(year_West_2_origcod), ncol = ncol(year_West_2_origcod))
  # g_{3, West}(year) is constant 0 for all years; already in original coding
  year_West_3 <- year_West_3_origcod
  # g_{1, East}(year)
  year_East_1_1991 <- - year_East_1_origcod[8,] + year_East_3_origcod[8,] +
    year_West_1_origcod[8,]
  year_East_1 <- t(apply(year_East_1_origcod, 1, "+", year_East_1_1991)) -
    year_East_3_origcod - year_West_1_origcod
  # g_{2, East}(year)
  year_East_2_1991 <- - year_East_2_origcod[8,] + year_East_3_origcod[8,] +
    year_West_2_origcod[8,]
  year_East_2 <- t(apply(year_East_2_origcod, 1, "+", year_East_2_1991)) -
    year_East_3_origcod - year_West_2_origcod
  # g_{3, East}(year) is constant 0 for all years
  year_East_3 <- matrix(0, nrow = nrow(year_East_3_origcod), ncol = ncol(year_East_3_origcod))


  ############# Add intercept to effects for plots on Bayes level ##############

  # Intercept + West_East
  West_East_intercept <- t(apply(West_East, 1, "+", intercept))

  # # Intercept + West_East + region
  # region_intercept <- region
  # for (i in 1:6) {
  #   West_East_ind <- ifelse(i < 5, 1, 2)
  #   region_intercept[i, ] <- West_East_intercept[West_East_ind, ] + region[i,]
  # }

  # Intercept + c_age
  c_age_intercept <- t(apply(c_age, 1, "+", intercept))

  # Intercept + West_East + c_age + West_East_c_age
  # West_East_intercept[1,] + c_age[1,] + West_East_c_age[1,] # West, 1
  # West_East_intercept[1,] + c_age[2,] + West_East_c_age[2,] # West, 2
  # West_East_intercept[1,] + c_age[3,] + West_East_c_age[3,] # West, 3
  # West_East_intercept[2,] + c_age[1,] + West_East_c_age[4,] # East, 1
  # ...
  West_East_c_age_intercept <- West_East_c_age
  for (i in 1:6) {
    West_East_ind <- ceiling(i / 3)
    cgroup_ind <- ifelse(i %% 3 > 0, i %% 3, 3)
    West_East_c_age_intercept[i, ] <- West_East_intercept[West_East_ind,] +
      c_age[cgroup_ind,] + West_East_c_age[i,]
  }

  # Intercept + year_intercept
  year_intercept_intercept <- t(apply(year_intercept, 1, "+", intercept))

  # Intercept + West_East + year_intercept + year_West_East
  year_West_intercept <- t(t(year_intercept_intercept) + West_East[1,] + t(year_West)) # Note that West_East[1,] = 0 = t(year_West)
  year_East_intercept <- t(t(year_intercept_intercept) + West_East[2,] + t(year_East))

  # Intercept + cgroup + year_intercept + year_c_age
  year_1_intercept <- t(t(year_intercept_intercept) + c_age[1,] + t(year_1))
  year_2_intercept <- t(t(year_intercept_intercept) + c_age[2,] + t(year_2))
  year_3_intercept <- t(t(year_intercept_intercept) + c_age[3,] + t(year_3)) # Note that c_age[3,] = 0 = t(year_3))

  # Intercept + West_East + c_age + West_East_c_age +
  # year_intercept + year_West_East + year_cgroup + year_West_East_cgroup
  year_West_1_intercept <- t(t(year_West_intercept) + c_age[1,] + West_East_c_age[1,] + t(year_1) + t(year_West_1)) # West, 1
  year_West_2_intercept <- t(t(year_West_intercept) + c_age[2,] + West_East_c_age[2,] + t(year_2) + t(year_West_2)) # West, 2
  year_West_3_intercept <- t(t(year_West_intercept) + c_age[3,] + West_East_c_age[3,] + t(year_3) + t(year_West_3)) # West, 3
  year_East_1_intercept <- t(t(year_East_intercept) + c_age[1,] + West_East_c_age[4,] + t(year_1) + t(year_East_1)) # East, 1
  year_East_2_intercept <- t(t(year_East_intercept) + c_age[2,] + West_East_c_age[5,] + t(year_2) + t(year_East_2)) # East, 2
  year_East_3_intercept <- t(t(year_East_intercept) + c_age[3,] + West_East_c_age[6,] + t(year_3) + t(year_East_3)) # East, 3

  ######################### Retransform to Bayes level #########################
  # Intercept
  intercept_pdf <- clr_inv(intercept, per_col = FALSE)
  intercept_d_pdf <- get_discrete_pdf(intercept_pdf)
  # Intercept + West_East
  West_East_intercept_pdf <- clr_inv(West_East_intercept, per_col = FALSE)
  West_East_intercept_d_pdf <- get_discrete_pdf(West_East_intercept_pdf)
  # # Intercept + West_East + region
  # region_intercept_pdf <- clr_inv(region_intercept, per_col = FALSE)
  # region_intercept_d_pdf <- get_discrete_pdf(region_intercept_pdf)
  # Intercept + c_age
  c_age_intercept_pdf <- clr_inv(c_age_intercept, per_col = FALSE)
  c_age_intercept_d_pdf <- get_discrete_pdf(c_age_intercept_pdf)
  # Intercept + West_East + c_age + West_East_c_age
  West_East_c_age_intercept_pdf <- clr_inv(West_East_c_age_intercept, per_col = FALSE)
  West_East_c_age_intercept_d_pdf <- get_discrete_pdf(West_East_c_age_intercept_pdf)
  # Intercept + year_intercept
  year_intercept_intercept_pdf <- clr_inv(year_intercept_intercept, per_col = FALSE)
  year_intercept_intercept_d_pdf <- get_discrete_pdf(year_intercept_intercept_pdf)
  # Intercept + West_East + year_intercept + year_West_East
  year_West_intercept_pdf <- clr_inv(year_West_intercept, per_col = FALSE)
  year_West_intercept_d_pdf <- get_discrete_pdf(year_West_intercept_pdf)
  year_East_intercept_pdf <- clr_inv(year_East_intercept, per_col = FALSE)
  year_East_intercept_d_pdf <- get_discrete_pdf(year_East_intercept_pdf)
  # Intercept + cgroup + year_intercept + year_cgroup
  year_1_intercept_pdf <- clr_inv(year_1_intercept, per_col = FALSE)
  year_1_intercept_d_pdf <- get_discrete_pdf(year_1_intercept_pdf)
  year_2_intercept_pdf <- clr_inv(year_2_intercept, per_col = FALSE)
  year_2_intercept_d_pdf <- get_discrete_pdf(year_2_intercept_pdf)
  year_3_intercept_pdf <- clr_inv(year_3_intercept, per_col = FALSE)
  year_3_intercept_d_pdf <- get_discrete_pdf(year_3_intercept_pdf)
  # Intercept + West_East + cgroup + West_East_c_age +
  # year_intercept + year_West_East + year_cgroup + year_West_East_cgroup
  year_West_1_intercept_pdf <- clr_inv(year_West_1_intercept, per_col = FALSE)
  year_West_1_intercept_d_pdf <- get_discrete_pdf(year_West_1_intercept_pdf)
  year_West_2_intercept_pdf <- clr_inv(year_West_2_intercept, per_col = FALSE)
  year_West_2_intercept_d_pdf <- get_discrete_pdf(year_West_2_intercept_pdf)
  year_West_3_intercept_pdf <- clr_inv(year_West_3_intercept, per_col = FALSE)
  year_West_3_intercept_d_pdf <- get_discrete_pdf(year_West_3_intercept_pdf)
  year_East_1_intercept_pdf <- clr_inv(year_East_1_intercept, per_col = FALSE)
  year_East_1_intercept_d_pdf <- get_discrete_pdf(year_East_1_intercept_pdf)
  year_East_2_intercept_pdf <- clr_inv(year_East_2_intercept, per_col = FALSE)
  year_East_2_intercept_d_pdf <- get_discrete_pdf(year_East_2_intercept_pdf)
  year_East_3_intercept_pdf <- clr_inv(year_East_3_intercept, per_col = FALSE)
  year_East_3_intercept_d_pdf <- get_discrete_pdf(year_East_3_intercept_pdf)


  ############### Childhood penalty (Diff in Diff) on clr-level ################
  # [f_{1, West} - f_{3, West}] - [f_{1, East} - f_{3, East}] # -> Each difference consists of cgroup + West_East_c_age + year_cgroup + year_West_East_cgroup
  childhood_penalty_West_1 <- year_West_1_intercept - year_West_3_intercept
  childhood_penalty_East_1 <- year_East_1_intercept - year_East_3_intercept
  diff_in_diff_1 <- childhood_penalty_West_1 - childhood_penalty_East_1
  # [f_{2, West} - f_{3, West}] - [f_{2, East} - f_{3, East}]
  childhood_penalty_West_2 <- year_West_2_intercept - year_West_3_intercept
  childhood_penalty_East_2 <- year_East_2_intercept - year_East_3_intercept
  diff_in_diff_2 <- childhood_penalty_West_2 - childhood_penalty_East_2

  # Transform to Bayes level
  diff_in_diff_1_pdf <- clr_inv(diff_in_diff_1, per_col = FALSE)
  diff_in_diff_1_d_pdf <- get_discrete_pdf(diff_in_diff_1_pdf)
  diff_in_diff_2_pdf <- clr_inv(diff_in_diff_2, per_col = FALSE)
  diff_in_diff_2_d_pdf <- get_discrete_pdf(diff_in_diff_2_pdf)

  ######### Return effects on clr level and + intercept on Bayes level #########

  return(list(G = G - 2, intercept = intercept,
              intercept_pdf = intercept_pdf, intercept_d_pdf = intercept_d_pdf,
              West_East = West_East, West_East_intercept_pdf = West_East_intercept_pdf,
              West_East_intercept_d_pdf = West_East_intercept_d_pdf,
              # region = region, region_intercept_pdf = region_intercept_pdf,
              # region_intercept_d_pdf = region_intercept_d_pdf,
              c_age = c_age, c_age_intercept_pdf = c_age_intercept_pdf,
              c_age_intercept_d_pdf = c_age_intercept_d_pdf,
              West_East_c_age = West_East_c_age,
              West_East_c_age_intercept_pdf = West_East_c_age_intercept_pdf,
              West_East_c_age_intercept_d_pdf = West_East_c_age_intercept_d_pdf,
              year_intercept = year_intercept,
              year_intercept_intercept_pdf = year_intercept_intercept_pdf,
              year_intercept_intercept_d_pdf = year_intercept_intercept_d_pdf,
              year_West = year_West, year_West_intercept_pdf = year_West_intercept_pdf,
              year_West_intercept_d_pdf = year_West_intercept_d_pdf,
              year_East = year_East, year_East_intercept_pdf  = year_East_intercept_pdf,
              year_East_intercept_d_pdf  = year_East_intercept_d_pdf,
              year_1 = year_1, year_1_intercept_pdf = year_1_intercept_pdf,
              year_1_intercept_d_pdf = year_1_intercept_d_pdf,
              year_2 = year_2, year_2_intercept_pdf = year_2_intercept_pdf,
              year_2_intercept_d_pdf = year_2_intercept_d_pdf,
              year_3 = year_3, year_3_intercept_pdf = year_3_intercept_pdf,
              year_3_intercept_d_pdf = year_3_intercept_d_pdf,
              year_West_1 = year_West_1, year_West_1_intercept_pdf = year_West_1_intercept_pdf,
              year_West_1_intercept_d_pdf = year_West_1_intercept_d_pdf,
              year_West_2 = year_West_2, year_West_2_intercept_pdf = year_West_2_intercept_pdf,
              year_West_2_intercept_d_pdf = year_West_2_intercept_d_pdf,
              year_West_3 = year_West_3, year_West_3_intercept_pdf = year_West_3_intercept_pdf,
              year_West_3_intercept_d_pdf = year_West_3_intercept_d_pdf,
              year_East_1 = year_East_1, year_East_1_intercept_pdf = year_East_1_intercept_pdf,
              year_East_1_intercept_d_pdf = year_East_1_intercept_d_pdf,
              year_East_2 = year_East_2, year_East_2_intercept_pdf = year_East_2_intercept_pdf,
              year_East_2_intercept_d_pdf = year_East_2_intercept_d_pdf,
              year_East_3 = year_East_3, year_East_3_intercept_pdf = year_East_3_intercept_pdf,
              year_East_3_intercept_d_pdf = year_East_3_intercept_d_pdf,
              diff_in_diff_1 = diff_in_diff_1, diff_in_diff_2 = diff_in_diff_2,
              diff_in_diff_1_pdf = diff_in_diff_1_pdf, diff_in_diff_2_pdf = diff_in_diff_2_pdf,
              diff_in_diff_1_d_pdf = diff_in_diff_1_d_pdf, diff_in_diff_2_d_pdf = diff_in_diff_2_d_pdf))
}

# # get_West_East_year_effects calculates predictions for West/East federal states over
# # the years. Arguments:
# # model_c: FDboost-object for the continuous model
# # model_d: FDboost-object for the discrete model
# # dat: underlying data for model estimation
# get_West_East_year_effects <- function(model_c, model_d, dat) {
#   ### Get estimated effects
#   West_East <- get_estimated_model_effects(model_c, model_d, dat, West_East = TRUE)$West_East
#
#   West <- matrix(rep(West_East[1, ], 33), nrow = 33, byrow = TRUE)
#   East <- matrix(rep(West_East[2, ], 33), nrow = 33, byrow = TRUE)
#
#   year <- 1984:2016
#   # Year West
#   newdat_West_c <- list(share = share, year = year,
#                        West_East = factor("West", levels = c("West", "East")))
#   newdat_West_d <- list(share_fac = c(0, 0.5, 1), year = year,
#                        West_East = factor("West", levels = c("West", "East")))
#   # year_intercept + year_West
#   year_West_c <- predict(model_c, newdata = newdat_West_c, which = 2) +
#     predict(model_c, newdata = newdat_West_c, which = 4)
#   year_West_d <- predict(model_d, newdata = newdat_West_d, which = 2) +
#     predict(model_d, newdata = newdat_West_d, which = 4)
#   year_West <- cbind(year_West_d[, 1], year_West_d[, 2] + year_West_c, year_West_d[, 3])
#   # Prediction per year for West Germany: intercept + West + year_intercept + year_West
#   year_West_pred <- West + year_West
#
#   # Year East
#   newdat_new_c <- list(share = share, year = year,
#                        West_East = factor("East", levels = c("West", "East")))
#   newdat_new_d <- list(share_fac = c(0, 0.5, 1), year = year,
#                        West_East = factor("East", levels = c("West", "East")))
#   # year_intercept + year_new
#   year_new_c <- predict(model_c, newdata = newdat_new_c, which = 2) + predict(model_c, newdata = newdat_new_c, which = 4)
#   year_new_d <- predict(model_d, newdata = newdat_new_d, which = 2) + predict(model_d, newdata = newdat_new_d, which = 4)
#   year_new <- cbind(year_new_d[, 1], year_new_d[, 2] + year_new_c, year_new_d[, 3])
#   # remove years before 1990
#   year_new[year < 1990,] <- NA
#   # Prediction per year for East Germany: intercept + East + year_intercept + year_East
#   year_East_pred <- East + year_East
#
#
#   ### Retransform to Bayes Hilbert space
#   year_West_pred_pdf <- exp(year_West_pred)
#   int_year_West_pred_pdf <- apply(year_West_pred_pdf, 1, integrate_num, estimation_points = share)
#   year_West_pred_pdf <- year_West_pred_pdf/int_year_West_pred_pdf
#
#   year_East_pred_pdf <- exp(year_East_pred)
#   int_year_East_pred_pdf <- apply(year_East_pred_pdf, 1, integrate_num, estimation_points = share)
#   year_East_pred_pdf <- year_East_pred_pdf/int_year_East_pred_pdf
#
#   return(list(share = share, year_West_pred_pdf = year_West_pred_pdf, year_East_pred_pdf = year_East_pred_pdf,
#               year_West_pred = year_West_pred, year_East_pred = year_East_pred))
# }
