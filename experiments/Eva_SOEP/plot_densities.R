### The following code creates all plots contained in Figure 4.1 (Estimated
### conditional densities)

################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/../.."))

def.par <- par(no.readonly = TRUE)

# Package to create nice color palettes
library(viridis)
# Load (amongst others) function to plot densities
source("./experiments/Eva_SOEP/help_functions.R")

path <- "./experiments/Eva_SOEP/Images/"
width <- 10

# Color gradient over the 33 observed years: from dark purple (1980s) to light
# purple (1990s), red/orange (2000s) to yellow (2016)
color <- plasma(33)
# Different line types for further differentiation
lty <- rep(c(1, 2, 4, 5), 9)

# color to indicate share level (0; (0, 1); 1)
color_breaks <- c("chartreuse2", "mediumblue", "red2")

# ################################################################################
# ################################# Plot Legends #################################
# ################################################################################
#
# # Year legend 2 columns
# pdf(paste0(path, "year_legend_v_2.pdf"), width = 2, height = 3.7)
# par(mar = c(0, 0, 0, 0) + 0.1)
# plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
# legend("topleft",xpd = TRUE, legend = 1984:2016, text.col = color, lty = lty,
#        cex = 1, ncol = 2, bty = "n", col = color, seg.len = 1.6,
#        title = "year", title.col = 1)
# dev.off()
#
# # # Year legend 1 column (If regions are included, legend for plot with all 6 regions)
# # pdf(paste0(path, "year_legend_v_1.pdf"), width = 0.95, height = 7)
# # par(mar = c(0, 0, 0.5, 0))
# # plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
# # legend("topleft",xpd = TRUE, legend = 1984:2016, text.col = color, lty = lty,
# #        cex = 1, ncol = 1, bty = "n", col = color, seg.len = 1.6,
# #        title = "year", title.col = 1)
# # dev.off()
#
# # Income share legend
# share_levels <- latex2exp::TeX(c("income share earned by the woman:",
#                                  "$s = 0$: non-working women",
#                                  "$s \\in (0, 1)$: dual-earner households",
#                                  "$s = 1$: single-earner women"))
# pdf(paste0(path, "legend_pred_hist_years_h_long.pdf"), width = 11.4, height = 0.5)
# par(mar = c(0, 0, 0, 0))
# plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
# text(x = c(0.09, 0.39, 0.66, 0.94), y = 0.5, share_levels, col = c(1, color_breaks))
# segments(x0 = c(0.27, 0.825), y0 = c(0.5, 0.5), x1 = c(0.28, 0.83), y1 = c(0.5, 0.5),
#          col = color_breaks[c(1, 3)])
# plotrix::draw.circle(0.525, 0.5, 0.0035, border = color_breaks[2])
# dev.off()


################################################################################
####################### Get data and functions for plots #######################
################################################################################

# Clr-predictions were calculated and retransformed in estimate_model.R
pred_clr_j <- readRDS("./experiments/Eva_SOEP/Objects/pred_clr_j.rds")
pred_j <- readRDS("./experiments/Eva_SOEP/Objects/pred_j.rds")
# apply(pred_clr_j, 2, integrate_num)
# apply(pred_j, 2, integrate_num)

# west_east <- c((33 * 3 + 1):(33 * 3 + 33* 3), (33 * 3 * 4 + 1):(33 * 3 * 4 + 33 * 3))
# ylim_west_east <- c(0, round(max(pred_j[, west_east], na.rm = TRUE), digits = 1))
ylim_all <- c(0, max(pred_j, na.rm = TRUE))
ylim_clr_all <- range(pred_clr_j, na.rm = TRUE)

# Get parameters for plots
get_param <- function(regions, up_down = 5) {
  n <<- length(regions) / 3
  c_age <<- rep(1:3, n)
  c_age_long <<- rep(c("0-6", "7-18", "other"), n)
  legend_text <<- paste0(#"region ",
    regions)
  legend_text_2 <<- paste0("c_age ", c_age_long)
  xaxt <<- rep(c(rep("n", n - 1),  "s"), each = 3)
  yaxt <<- rep(c("s", "n", "n"), n)
  mar_inner <<- rep(list(c(0, 4, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 4)), n - 2)
  mar <<- c(list(c(0, 4, up_down, 0), c(0, 0, up_down, 0), c(0, 0, up_down, 4)),
            mar_inner,
            list(c(up_down, 4, 0, 0), c(up_down, 0, 0, 0), c(up_down, 0, 0, 4)))
  all_regions <- c("West", "East") # c("northwest", "west", "southwest", "south", "east", "northeast")
  regions_pos <- unique(which(all_regions %in% regions))
  j <<- rep(regions_pos, each = 3)
  y_labels <- c(seq(0, 0.8, by = 0.2), "")
  y_labels <<- rep(list(y_labels, NULL, NULL), n)
}

# Function to plot several densities arranged in a matrix, using plot_densities
# contained in ./Functions/help_functions.R
plot_dens_i <- function (i, dens_matrix, ylab, ylim, cex, cex_legend,
                         x_legend = "topright", y_legend = NULL,
                         y_at = NULL, y_labels = NULL, ...) {
  plot_densities(dens_matrix, path = path, ylim = ylim, pdf = FALSE, main = "",
                 regions = regions[i], c_age = c_age[i], mar = mar[[i]],
                 xaxt = xaxt[i], yaxt = yaxt[i], las = 1, ylab = ylab, cex = cex,
                 xlab = "", y_at = y_at, y_labels = y_labels, ...)
  legend(x = x_legend, y = y_legend, legend = c(legend_text[i], legend_text_2[i]), bty = "n", cex = cex_legend)
}

# Function to get corresponding column-numbers (years) in the density/prediction
# matrix for given integers corresponding to region (j) and c_age (k)
get_year_vec <- function(j, k, all_years = TRUE) {
  intercept_West_East <- ifelse(all_years, 1, ifelse(j < 5,  1, 8))
  first_year <- (j - 1) * 3 * 33 + (k - 1) * 33 + intercept_West_East
  last_year <- (j - 1) * 3 * 33 + k * 33
  return(first_year:last_year)
}
# Function to plot relative frequencies, using plot_histograms_per_year
# contained in ./Functions/help_functions.R
plot_j_k <- function(dens_matrix, j, k, mar = NULL, xaxt = "s", yaxt = "s",
                     regions = c("northwest", "west", "southwest", "south", "east", "northeast"),
                     y_at = NULL, y_labels = NULL, plot_type = "p", cex_legend, c_age,
                     x_legend = 1999, y_legend = 1.15, ...) {
  initial_j_k <- apply(dens_matrix[, get_year_vec(j, k)], 2, integrate_num, partial = TRUE)
  plot_histograms_per_year(hist_data = initial_j_k, pdf = FALSE, col = color_breaks,
                           ylim = c(0, 1), y_at = y_at, y_labels = y_labels,
                           text_labels = NULL, las = 1, plot_type = plot_type,
                           ylab = "", mar = mar, xaxt = xaxt, yaxt = yaxt, ...)
  legend(x = x_legend, y = y_legend, bty = "n", cex = cex_legend, xjust = 0.5,
         legend = paste0(regions[j], ", ", c_age_long[k]))
}

################################################################################
########### Plot estimated (clr-)densities for West and East Germany ###########
################################################################################

regions <- rep(c("West", "East"), each = 3)
get_param(regions, up_down = 4.2)
height <- 4

### Predictions on clr-level
pdf(paste0(path, "clr_pred_West_East.pdf"), width = width, height = height)
layout(matrix(seq_along(regions), ncol = 3, byrow = TRUE), widths = c(1, 0.85, 1))
sapply(seq_along(regions), plot_dens_i, dens_matrix = pred_clr_j, ylab = "",
       ylim = ylim_clr_all, cex = 0.8, cex_legend = 1, x_legend = "bottom")
mtext(text = "predicted densities", side = 2, line = 39.5, at = ylim_clr_all[2])
mtext(text = "s: income share earned by the woman", side = 1, line = 3, at = -0.6)
dev.off()

### Predictions on Bayes level (Figure 4.1)
# Densities (upper part of Figure 4.1)
pdf(paste0(path, "pred_West_East.pdf"), width = width, height = height)
layout(matrix(seq_along(regions), ncol = 3, byrow = TRUE), widths = c(1, 0.85, 1))
sapply(seq_along(regions), plot_dens_i, dens_matrix = pred_j, ylab = "",
       ylim = ylim_all, cex = 0.8, cex_legend = 1, x_legend = 0.53, y_legend = 3,
       y_at = seq(0, 2.5, by = 0.5))
mtext(text = "predicted densities", side = 2, line = 39.5, at = ylim_all[2])
mtext(text = "s: income share earned by the woman", side = 1, line = 3, at = -0.6)
dev.off()

# Relative frequencies per year (lower part of Figure 4.1)
height <- 2.8
pdf(paste0(path, "pred_West_East_hist_years.pdf"), width = width, height = height)
layout(matrix(seq_along(regions), ncol = 3, byrow = TRUE), widths = c(1, 0.85, 1))
sapply(seq_along(regions), function(i) plot_j_k(pred_j, j[i], c_age[i], mar = mar[[i]],
                                                plot_type = "l", lty = c(1, 3, 2), xaxt = xaxt[i],
                                                yaxt = yaxt[i], y_labels = y_labels[[i]],
                                                cex = 0.8, cex_legend = 1, xlab = "",
                                                y_legend = 1.15, regions = unique(regions)))
mtext(text = "predicted relative frequency", side = 2, line = 39.5, at = 1) # 3.2
mtext(text = "year", side = 1, line = 3, at = 1982 - length(1983:2016) / 2)
dev.off()

par(def.par)

################################################################################
######## Plot underlying histograms with estimated conditional density #########
########                    (not included in paper)                    #########
################################################################################

# For each covariate combination: Compute histogram on relative (i.e., density)
# scale and add estimated density
hist_data <- readRDS("./experiments/Eva_SOEP/Objects/hist_data.rds")

share <- hist_data$share
year <- c(rep(1984:2016, 3), rep(1991:2016, 3))
c_age <- c(rep(c("0-6", "7-18", "other"), each = length(1984:2016)),
           rep(c("0-6", "7-18", "other"), each = length(1991:2016)))
West_East <- c(rep("West", 3 * length(1984:2016)), rep("East", 3 * length(1991:2016)))
color <- c(rep(color, 3), rep(color[-(1:7)], 3))

# Remove artificial NA-columns (corresponding to 1984-1990 in East Germany) for
# single plots
f_hat <- pred_j[, -(rep(c(100, 133, 166), each = 7) + 0:6)] # which(apply(pred_j, 2, function(x) any(is.na(x))))

color_hist <- rep("grey", 100)
color_hist[51] <- "red" # bar containing share of 0.5 is highlighted in red
                        # (Sonja told me to check, because in histogram of complete
                        # data, there is a spike at 0.5)
library(manipulate)
manipulate({
  barplot(hist_data$dual_earner_dens[,i], width = max(diff(share)), space = 0,
          main = paste0(West_East[i], " ", year[i], ", c_age: ", c_age[i]),
          ylim = c(0, max(hist_data$zeros_dens[i], hist_data$dual_earner_dens[,i],
                          hist_data$ones_dens[i], f_hat[, i], na.rm = TRUE) + 0.010),
          col = color_hist)
  axis(1)
  points(c(0, 1), c(hist_data$zeros_dens[i], hist_data$ones_dens[i]))
  if (!any(is.na(f_hat[, i]))) {
    lines(share[2:101], f_hat[2:101, i], type = "l", col = color[i])
    points(c(0, 1), f_hat[c(1, 102), i], col = color[i], pch = 4)
  }
}, i = slider(1, length(year)))

# Estimated conditional densities look reasonable
# Regarding 0.5, there are several covariate combinations with a high bar value at
# 0.5, but (mostly) not out of scale compared to other bars in the same histograms.
# On the other hand, there are also histograms with a low bar value at 0.5 and/or
# high bar values at other share values.
# In my opinion, this shows that there is no special treatment of 0.5 necessary.

pdf(paste0(path, "histograms_West.pdf"), width = width, height = 75)
layout(matrix(1:99, nrow = 33, ncol = 3))
for(i in 1:99) {
  barplot(hist_data$dual_earner_dens[,i], width = max(diff(share)), space = 0,
          main = paste0(West_East[i], " ", year[i], ", c_age: ", c_age[i]),
          ylim = c(0, max(hist_data$zeros_dens[i], hist_data$dual_earner_dens[,i],
                          hist_data$ones_dens[i], f_hat[, i], na.rm = TRUE) + 0.010),
          col = color_hist)
  axis(1)
  points(c(0, 1), c(hist_data$zeros_dens[i], hist_data$ones_dens[i]))
  if (!any(is.na(f_hat[, i]))) {
    lines(share[2:101], f_hat[2:101, i], type = "l", col = color[i])
    points(c(0, 1), f_hat[c(1, 102), i], col = color[i], pch = 4)
  }
}
dev.off()

pdf(paste0(path, "histograms_East.pdf"), width = width, height = 60)
layout(matrix(1:78, nrow = 26, ncol = 3))
for(i in 100:177) {
  barplot(hist_data$dual_earner_dens[,i], width = max(diff(share)), space = 0,
          main = paste0(West_East[i], " ", year[i], ", c_age: ", c_age[i]),
          ylim = c(0, max(hist_data$zeros_dens[i], hist_data$dual_earner_dens[,i],
                          hist_data$ones_dens[i], f_hat[, i], na.rm = TRUE) + 0.010),
          col = color_hist)
  axis(1)
  points(c(0, 1), c(hist_data$zeros_dens[i], hist_data$ones_dens[i]))
  if (!any(is.na(f_hat[, i]))) {
    lines(share[2:101], f_hat[2:101, i], type = "l", col = color[i])
    points(c(0, 1), f_hat[c(1, 102), i], col = color[i], pch = 4)
  }
}
dev.off()
