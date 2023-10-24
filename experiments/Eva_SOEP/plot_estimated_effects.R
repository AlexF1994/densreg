### The following code creates all plots (except for region, which is so far not
### included as covariate in the model) corresponding to the plots contained in
### Figures 2-5 (as well as Figures E.9-E.16 in the appendix) in the first paper.

################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/../.."))

def.par <- par(no.readonly = TRUE)

# Package to create nice color palettes
library(RColorBrewer)
# Package to convert latex strings to R expressions
library(latex2exp)
# Load (amongst others) functions to compute and plot estimated effects
source("./experiments/Eva_SOEP/help_functions.R")
# load self-written smoother for mixed reference measure
source("./R/mixed_density_smooth.R")

path <- "./experiments/Eva_SOEP/Images/"

model_soep <- readRDS("./experiments/Eva_SOEP/Objects/model_all_covariates_but_region.rds")

### Confidence region?
# library(stringr)
# effect_names <- unique(str_extract(names(model_soep$coefficients), regex(".*(?=\\.[0-9]+)")))
# effect_names <- effect_names[-which(is.na(effect_names))]
#
# intercept_ind <- which(grepl(paste0(effect_names[1], "."), names(model_soep$coefficients), fixed = TRUE))
#
# # X <- model.matrix(model)[, intercept_ind]
# theta_hat <- model$coefficients[intercept_ind]
#
# # f_hat_clr <- c(X %*% theta_hat)
# # f_hat <- exp(f_hat_clr) / mean(exp(f_hat_clr))
#
# # theta ~ N(theta_hat, Vc)
# # => CR = {theta : (theta - theta_hat)^T V^{-1} (theta - theta_hat) <= chi^2_{1-alpha}(K_T)}
# # An effect is significant, if 0 \notin CR
#
# # Das entspricht noch nicht unserer Theorie! In Simulation (woraus der Code
# # kopiert ist) gab es keine Kovariablen, deshalb war theta = theta_T; Hier ist
# # aber vermutlich die Transformation wie im Paper in 2.4 beschreiben, nötig
# # -> Noch mal überprüfen/nachvollziehen, wenn wach und geistig fit
# check_significance <- function(model, indices, type = c("Vc", "Vp"), alpha = 0.05) {
#   type <- match.arg(type)
#   if (type == "Vc") {
#     V <- model$Vc[indices, indices, drop = FALSE]
#   } else {
#     V <- model$Vp[indices, indices, drop = FALSE]
#   }
#   theta_hat <- model$coefficients[indices]
#   theta_diff <- matrix(0 - theta_hat, ncol = 1)
#   V_inv <- solve(V)
#   chi_statistic_V <- t(theta_diff) %*% V_inv %*% theta_diff
#   check_coverage_V <- as.numeric(chi_statistic_V) <= qchisq(1-alpha, df = length(indices)) # is 0 covered?
#
#   return(list(significance = !(check_coverage_V), statistic = chi_statistic_V,
#               quantile = qchisq(1-alpha, df = length(indices))))
# }
#
# library(stringr)
# effect_names <- sapply(model_soep$smooth, function(s) s$label) #unique(str_extract(names(model_soep$coefficients), regex(".*(?=\\.[0-9]+)")))
# # effect_names <- effect_names[-which(is.na(effect_names))]
#
# # i = 4 (ti(share):c_age1): Covariance not invertible (system is computationally singular: reciprocal condition number = 1.5754e-17)
# # i = 5 (ti(share):West_East_c_ageWest_2): Covariance not invertible (system is exactly singular: U[12,12] = 0)
# # i = 6 (ti(share):West_East_c_ageWest_1): Covariance not invertible (system is exactly singular: U[12,12] = 0)
# # i = 8 (ti(share):West_East_c_ageEast_2): Covariance not invertible (system is exactly singular: U[12,12] = 0)
# # i = 9 (ti(share):West_East_c_ageEast_1): Covariance not invertible (system is exactly singular: U[12,12] = 0)
# significance <- lapply(seq_along(effect_names)[-c(4:6, 8:9)], function(i) check_significance(model_soep, which(grepl(paste0(effect_names[i], "."), names(model_soep$coefficients), fixed = TRUE))))
# names(significance) <- paste0(seq_along(effect_names)[-c(4:6, 8:9)], ": ", effect_names[-c(4:6, 8:9)])
# # Not significant:
# # i = 11 (ti(share,syear):West_EastEast)
# # i = 12 (ti(share,syear):c_age2)
# # i = 13 (ti(share,syear):c_age1)
# # i = 14 (ti(share,syear):West_East_c_ageWest_2)
# # i = 15 (ti(share,syear):West_East_c_ageWest_1)
# # i = 16 (ti(share,syear):West_East_c_ageEast_3)
# # i = 17 (ti(share,syear):West_East_c_ageEast_2)
# # i = 18 (ti(share,syear):West_East_c_ageEast_1)

# Compute effects to plot
effects <- get_estimated_model_effects(model_soep, G = 100)
saveRDS(effects, "./experiments/Eva_SOEP/Objects/estimated_effects.rds")

# Color gradient over the 33 observed years: from dark purple (1980s) to light
# purple (1990s), red/orange (2000s) to yellow (2016)
color <- viridis::plasma(33)
# Different line types for further differentiation
lty_year <- rep(c(1, 2, 4, 5), 9)

################################################################################
################### Plots for Section 4 (main part of paper) ###################
################################################################################

############ Figure 2 (and Figure E.11 in appendix): c_age effects #############
c_age_col <- brewer.pal(n = 3, "Dark2")
c_age_legend <- c("c_age = 0-6", "c_age = 7-18", "c_age = other")

# Plot c_age on clr-level
plot_effects(effects$c_age, col = c_age_col, path = path, abline_h = NULL,
             ylab = expression(paste(clr, " [", hat(beta)["c_age"], "]")), las = 1,
             legend = c_age_legend, legend_x = "bottomleft", legend_cex = 0.7)

# Plot intercept \oplus c_age on Bayes level
ylab_intercept_cgroup <- expression(paste("", "", hat(paste("", beta)),
                                           phantom()[{paste("0")}], "", paste(" "), "",
                                           symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                           phantom()[{paste("", "c_age")}], ""))

plot_effects(effects$c_age_intercept_pdf, col = c_age_col, path = path,
             effect_matrix_d = effects$c_age_intercept_d_pdf, pch = 1,
             ylab = ylab_intercept_cgroup, las = 1, legend_cex = 0.7,
             legend = c_age_legend, legend_x = "topright")

#### Figure 3/4: Expected densities for West/East Germany (Figure 3/4) for #####
####       1984, 1991, 2003, and 2016 for the different c_age values       #####

# Get Intercept \oplus year \oplus year_West_East \oplus year_c_age
# \oplus year_West_East_c_age on Bayes level
year_West_c_age_intercept <- list(effects$year_West_3_intercept_pdf,
                                 effects$year_West_1_intercept_pdf,
                                 effects$year_West_2_intercept_pdf)
year_East_c_age_intercept <- list(effects$year_East_3_intercept_pdf,
                                 effects$year_East_1_intercept_pdf,
                                 effects$year_East_2_intercept_pdf)
year_West_c_age_intercept_d <- list(effects$year_West_3_intercept_d_pdf,
                                   effects$year_West_1_intercept_d_pdf,
                                   effects$year_West_2_intercept_d_pdf)
year_East_c_age_intercept_d <- list(effects$year_East_3_intercept_d_pdf,
                                   effects$year_East_1_intercept_d_pdf,
                                   effects$year_East_2_intercept_d_pdf)

ylim <- c(0, round(max(sapply(year_West_c_age_intercept, max, na.rm = TRUE),
                       sapply(year_East_c_age_intercept, max, na.rm = TRUE)), 1))
yaxt <- c("s", "n", "n")
xlab <- c("", "s: income share earned by the woman", "")
legend_text <- paste0("c_age = ", c("other", "0-6", "7-18"))
mar <- c(list(c(5, 11, 5, 0), c(5, 0, 5, 0), c(5, 0, 5, 11)))

# Labels for Figure 3
y_labels_West_1 <- expression(paste(hat(f)[{paste("", "c_age, West, year")}], " := "))
y_labels_West_2 <- expression(paste(hat(paste("", beta)), phantom()[{paste("0")}], "", paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "West")}], paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age")}], paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age, West")}], paste(" "), ""))
y_labels_West_3 <- expression(paste(symbol("\305"), paste(" "), "", hat(paste("", g)),  "", (year), paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "West")}], "", (year), paste(" "), ""))
y_labels_West_4 <- expression(paste(symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age")}], "", (year), paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age, West")}], "", (year)))

# Labels for Figure 4
y_labels_East_1 <- expression(paste(hat(f)[{paste("", "c_age, East, year")}], " := "))
y_labels_East_2 <- expression(paste("", "", hat(paste("", beta)), phantom()[{paste("0")}], "", paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "East")}], paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age")}], paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age, East")}], paste(" "), ""))
y_labels_East_3 <- expression(paste(symbol("\305"), paste(" "), "", hat(paste("", g)), "", (year), paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "East")}], "", (year), paste(" "), ""))
y_labels_East_4 <- expression(paste(symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age")}], "", (year), paste(" "), "",
                                   symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age, East")}], "", (year)))
# Position of Labels
at <- 0.9

# Only consider selected years: 1984, 1991, 2003, 2016
selected_years <- c(1984, 1991, 2003, 2016) - 1983
year_West_c_age_intercept_sel <- lapply(year_West_c_age_intercept, function(X) X[selected_years,])
year_East_c_age_intercept_sel <- lapply(year_East_c_age_intercept, function(X) X[selected_years,])
year_West_c_age_intercept_d_sel <- lapply(year_West_c_age_intercept_d, function(X) X[selected_years,])
year_East_c_age_intercept_d_sel <- lapply(year_East_c_age_intercept_d, function(X) X[selected_years,])

# Figure 3
pdf(paste0(path, "estimated_sel_year_West_c_age_intercept_h.pdf"), width = 10.5, height = 3)
layout(matrix(1:3, ncol = 3, byrow = TRUE), widths = c(1, 0.64, 1))
sapply(1:3, function(i) #par(mar = mar[[i]])
  plot_effects(year_West_c_age_intercept_sel[[i]], col = color[selected_years],
               lty = lty_year[selected_years], pdf = FALSE,
               ylim = ylim, mar = mar[[i]], # lwd = 1,
               effect_matrix_d = year_West_c_age_intercept_d_sel[[i]], pch = 1,
               legend = legend_text[[i]], legend_x = "topright",
               legend_col = 1, legend_cex = 1.2,
               # main = "Est. Effect of year for c_age 1 (clr-level)",
               ylab = "", cex.axis = 1, cex.lab = 1.2,
               xlab = xlab[i], #xaxt = "n",
               yaxt = yaxt[i], las = 1, abline_h = NULL))
mtext(text = y_labels_West_1, side = 2, line = 46.8, at = at, cex = 0.9)
mtext(text = y_labels_West_2, side = 2, line = 44.8, at = at, cex = 0.9)
mtext(text = y_labels_West_3, side = 2, line = 42.8, at = at, cex = 0.9)
mtext(text = y_labels_West_4, side = 2, line = 40.8, at = at, cex = 0.9)
legend(x = 1.1, y = 1.5, legend = selected_years + 1983, col = color[selected_years],
       text.col = color[selected_years], lty = lty_year[selected_years],
       title = "year", title.col = 1, xpd = TRUE, bty = "n", cex = 1.2)
dev.off()

# Figure 4
pdf(paste0(path, "estimated_sel_year_East_c_age_intercept_h.pdf"), width = 10.5, height = 3)
layout(matrix(1:3, ncol = 3, byrow = TRUE), widths = c(1, 0.64, 1))
sapply(1:3, function(i) #par(mar = mar[[i]])
  plot_effects(year_East_c_age_intercept_sel[[i]], col = color[selected_years],
               lty = lty_year[selected_years], pdf = FALSE,
               ylim = ylim, mar = mar[[i]], # lwd = 1,
               effect_matrix_d = year_East_c_age_intercept_d_sel[[i]], pch = 1,
               legend = legend_text[[i]], legend_x = "topright",
               legend_col = 1, legend_cex = 1.2,
               # main = "Est. Effect of year for c_age 1 (clr-level)",
               ylab = "", cex.axis = 1, cex.lab = 1.2,
               xlab = xlab[i], #xaxt = "n",
               yaxt = yaxt[i], las = 1, abline_h = NULL))
mtext(text = y_labels_East_1, side = 2, line = 46.8, at = at, cex = 0.9)
mtext(text = y_labels_East_2, side = 2, line = 44.8, at = at, cex = 0.9)
mtext(text = y_labels_East_3, side = 2, line = 42.8, at = at, cex = 0.9)
mtext(text = y_labels_East_4, side = 2, line = 40.8, at = at, cex = 0.9)
legend(x = 1.1, y = 1.5, legend = selected_years[-1] + 1983, col = color[selected_years[-1]],
       text.col = color[selected_years[-1]], lty = lty_year[selected_years[-1]],
       title = "year", title.col = 1, xpd = TRUE, bty = "n", cex = 1.2)
dev.off()

#### Figure 5: Log Odds of the West-East gap in the childhood penalty (DiD) ####
width <- 7
height <- 4
mar_up_down <- 4.2
mar_between <- 0.4
mar_left_right <- 5
x_axis <- rep(c(FALSE, TRUE), each = 2)
y_axis <- rep(c(TRUE, FALSE), 2)

mar <- c(list(c(mar_between, mar_left_right, mar_up_down, mar_between),
              c(mar_between, mar_between, mar_up_down, mar_left_right)),
         list(c(mar_up_down, mar_left_right, mar_between, mar_between),
              c(mar_up_down, mar_between, mar_between, mar_left_right)))
main <- c("1991", "2016", "", "")
legend <- list(NULL, "c_age\n= 0-6", NULL, "c_age\n= 7-18")

diff_in_diff_odds_all <- list(effects$diff_in_diff_1[selected_years[2],],
                              effects$diff_in_diff_1[selected_years[4],],
                              effects$diff_in_diff_2[selected_years[2],],
                              effects$diff_in_diff_2[selected_years[4],])
diff_in_diff_odds_all_d <- lapply(diff_in_diff_odds_all, function(x) integrate_num(x, partial = TRUE))

diff_in_diff_odds_all <- lapply(diff_in_diff_odds_all,
                                function (D) sapply(D, function(i) i - D))
diff_in_diff_odds_all_d <- lapply(diff_in_diff_odds_all_d,
                                  function (D) sapply(D, function(i) i - D))

collim_complete <- max(abs(unlist(diff_in_diff_odds_all)))
collim_complete <- ceiling(collim_complete * 10) / 10
collim <- 1.5
# colorRampPalette interpolates color palettes, rev inverts the palette
n_image_col <- 2 * collim * 10 + 2
image_col <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(n_image_col))
image_breaks <- c(-2, seq(-collim, collim, length = n_image_col - 1), 2)
image_breaks_contour <- seq(-collim_complete, collim_complete, length = (2 * collim_complete * 10) + 1)
# image_breaks <- seq(-2, 2, by = 0.2)

pdf(paste0(path, "estimated_cp_diff_in_diff_odds.pdf"), width = width, height = height)
layout(matrix(1:4, ncol = 2, byrow = TRUE))
sapply(seq_along(diff_in_diff_odds_all), function(i)
  plot_odds(diff_in_diff_odds_all[[i]], diff_in_diff_odds_all_d[[i]],
            image_col = image_col,
            image_breaks = image_breaks, mar = mar[[i]], contour_cex = 0.6,
            x_axis = x_axis[i], y_axis = y_axis[i],
            legend = legend[[i]], legend_year = main[i], n_border_left_right = 2,
            n_border_up_down = 4, legend_cex = 0.8, axis_cex = 0.8))
mtext(text = "s: income share earned by the woman", side = 2,
      line = 20, at = 1.05,
      cex = 0.75)
mtext(text = "t: income share earned by the woman", side = 1,
      line = 3, at = -0.05,
      cex = 0.75)
dev.off()

# legend for heatmap vertical
# image_col_mod <- c(rep(image_col[1], 4), image_col,
#                    rep(image_col[length(image_col)], 4))
image_col_mod <- image_col
pdf(paste0(path, "cp_diff_in_diff_legend_v.pdf"), width = 0.95, height = 4)
par(mar = c(0, 2.5, 0, 0))
barplot(height = rep(1, length(image_col_mod)), width = 1, space = 0, border = NA,
        col = image_col_mod, horiz = TRUE, xaxt = "n")
legend_at <- seq(0, length(image_breaks) - 1, by = 5)
axis(2, at = legend_at, labels = image_breaks[legend_at + 1], las = 1,
     lwd = 0, lwd.ticks = 1, pos = 0)
dev.off()

################################################################################
################## Plot remaining estimated effects (Appendix) #################
################################################################################

########### Plot West_East (East_West) effects: Figure E.9 in appendix ###########

West_East_col <- brewer.pal(3, "Set1")[2:1]

# Plot West/East on clr-level
# Note that the function plot_effects saves the plot as pdf (under the given
# path) per default; If that's not desired, set pdf = FALSE
plot_effects(effects$West_East, col = West_East_col, path = path, abline_h = NULL,
             ylab = expression(paste("clr [",hat(beta)[West_East], "]")),
             legend = c("West_East = West", "West_East = East"), las = 1,
             legend_x = "bottomright", legend_cex = 0.7)

# Plot intercept \oplus West/East on Bayes level

# symbol("\305") corresponds to \oplus
ylab_intercept_West_East <- expression(paste("", "", hat(paste("", beta)),
                                           phantom()[{paste("0")}], "", paste(" "), "",
                                           symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                           phantom()[{paste("", "West_East")}], ""))

plot_effects(effects$West_East_intercept_pdf, col = West_East_col, path = path,
             pch = 1, effect_matrix_d = effects$West_East_intercept_d_pdf,
             ylab = ylab_intercept_West_East, las = 1, legend_cex = 0.7,
             legend = c("West_East = West", "West_East = East"), legend_x = "topright")

# ################# Plot region effects (Figure E.10 in appendix) ################
#
# region_col <- brewer.pal(9, "Set1")[c(3, 2, 1, 4, 5, 9)]
# regions <- c("northwest", "west", "southwest", "south", "east", "northeast")
# region_legend <- paste0("region = ", regions, " (", c(rep("West", 4), rep("East", 2)), ")")
#
# # Plot region on clr-level
# plot_effects(effects$region, col = region_col, path = path,
#              ylab = expression(paste(clr, " [", hat(beta)[region], "]")), las = 1,
#              legend = region_legend, legend_x = "bottomleft", legend_cex = 0.7,
#              legend_ncol = 2)
#
# # Plot intercept \oplus West_East \oplus c_age on Bayes level
# ylab_intercept_region <- expression(paste("", "", hat(paste("", beta)),
#                                           phantom()[{paste("0")}], "", paste(" "), "",
#                                           symbol("\305"), paste(" "), "", hat(paste("", beta)),
#                                           phantom()[{paste("", "West_East")}], "", paste(" "), "",
#                                           symbol("\305"), paste(" "), "", hat(paste("", beta)),
#                                           phantom()[{paste("", "region")}], ""))
# plot_effects(effects$region_intercept_pdf, col = region_col, path = path,
#              effect_matrix_d = effects$region_intercept_d_pdf, pch = 1,
#              ylab = ylab_intercept_region, las = 1, legend_cex = 0.7,
#              legend = region_legend, legend_x = "topright")

### Plot interaction effects for West_East and c_age (Figure E.12 in appendix) ###

interaction_col <- c(brewer.pal(n = 3, "Dark2"), brewer.pal(n = 3, "Accent")[c(1, 3, 2)])
interaction_legend <- paste0("c_age = ", c("0-6", "7-18", "other"),
                             ", West_East = ", c(rep("West", 3), rep("East", 3)))

# Plot interaction on clr-level
plot_effects(effects$West_East_c_age, col = interaction_col, path = path,
             ylab = expression(paste(clr, " [", hat(beta)["c_age, West_East"], "]")),
             legend = interaction_legend, legend_x = "bottom", legend_ncol = 2,
             legend_cex = 0.7, abline_h = NULL, las = 1)

# Plot intercept \oplus West_East \oplus c_age \oplus West_East_c_age on Bayes level
ylab_West_East_c_age_1 <- expression(paste("", "", hat(paste("", beta)),
                                          phantom()[{paste("0")}], "", paste(" "), "",
                                          symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                          phantom()[{paste("", "West_East")}], "", paste(" "), "",
                                          symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                          phantom()[{paste("", "c_age")}], ""))
ylab_West_East_c_age_2 <- expression(paste(symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                          phantom()[{paste("", "c_age, West_East")}], ""))
at_West_East_c_age <- max(effects$West_East_c_age_intercept_pdf, na.rm = TRUE)/2

pdf(paste0(path, "estimated_West_East_c_age_intercept.pdf"), width = 5, height = 3)
plot_effects(effects$West_East_c_age_intercept_pdf, col = interaction_col,
             path = path, pdf = FALSE,  pch = 1, ylab = "",
             effect_matrix_d = effects$West_East_c_age_intercept_d_pdf, las = 1,
             legend_cex = 0.7, legend = interaction_legend, legend_x = "topright")
mtext(text = ylab_West_East_c_age_1, side = 2, line = 3.8, at = at_West_East_c_age)
mtext(text = ylab_West_East_c_age_2, side = 2, line = 2.5, at = at_West_East_c_age)
dev.off()

################# Plot year effects (Figure E.13 in appendix) ##################

# Legend for years horizontal
pdf(paste0(path, "year_legend.pdf"), width = 15.3, height = 1)
par(mar = c(1.3,0,0,0) + 0.1, mgp = c(0, 0.2,0), cex = 1.4)
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("topleft",xpd = TRUE, legend = 1984:2016, text.col = color, lty = lty_year,
       cex = 1, ncol = 11, bty = "n", col = color)
dev.off()

# Plot year on clr-level
plot_effects(effects$year_intercept, col = color, lty = lty_year, path = path,
             abline_h = NULL, las = 1,
             ylab = expression(paste(clr, " [", hat(g)(year), "]")))

# Plot intercept \oplus year on Bayes level
plot_effects(effects$year_intercept_intercept_pdf, col = color,  lty = lty_year,
             path = path, effect_matrix_d = effects$year_intercept_intercept_d_pdf,
             pch = 1, ylab = TeX("$\\hat{\\beta}_0\\, \\oplus\\, \\hat{g}(year)$"),
             las = 1)

############ Plot year effect for West_East (Figure E.14 in appendix) ############

# plot on clr-level arranged as matrix vertical
le <- 5.3
mar <- c(list(c(0, le, 5, 1), c(5, le, 0, 1)))
y_labels <- expression(paste(clr, " [", hat(g)[West_East](year), "]"))
xaxt <- c("n", "s")
legend_text <- paste0("West_East = ", c("West", "East"))
year_West_East <- list(effects$year_West, effects$year_East)
year_West_East_d <- list(effects$year_West_d, effects$year_East_d)
ylim_clr_year_West_East <- range(c(effects$year_West, effects$year_East), na.rm= TRUE)

pdf(paste0(path, "estimated_year_West_East_clr.pdf"), width = 5, height = 5.7)
layout(matrix(1:2, ncol = 1))
sapply(1:2, function(i) plot_effects(year_West_East[[i]], col = color, lty = lty_year,
                                     pdf = FALSE, ylim = ylim_clr_year_West_East,
                                     mar = mar[[i]], legend = legend_text[[i]],
                                     legend_x = "top", legend_col = 1, ylab = "",
                                     xaxt = xaxt[i], las = 1, abline_h = NULL))
mtext(text = y_labels, side = 2, line = 3, at = max(year_West_East[[2]], na.rm = TRUE))
dev.off()

# Plot Intercept \oplus West_East \oplus year \oplus year_West_East on Bayes level
y_labels <- expression(paste("", "", hat(paste("", beta)), phantom()[{paste("0")}], "", paste(" "), "",
                             symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "West_East")}], paste(" "), "",
                             symbol("\305"), paste(" "), "", hat(paste("", g)),  "", (year), paste(" "), "",
                             symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "West_East")}], "", (year)))
xaxt <- c("n", "s")
legend_text <- paste0("West_East = ", c("West", "East"))
year_West_East_intercept <- list(effects$year_West_intercept_pdf, effects$year_East_intercept_pdf)
year_West_East_intercept_d <- list(effects$year_West_intercept_d_pdf, effects$year_East_intercept_d_pdf)
ylim <- c(0, round(max(unlist(lapply(year_West_East_intercept, max, na.rm = TRUE))), 1))

pdf(paste0(path, "estimated_year_West_East_intercept.pdf"), width = 5, height = 5.7)
layout(matrix(1:2, ncol = 1))
sapply(1:2, function(i) plot_effects(year_West_East_intercept[[i]], col = color,
                                     lty = lty_year, pdf = FALSE, ylim = ylim,
                                     mar = mar[[i]], pch = 1, ylab = "",
                                     effect_matrix_d = year_West_East_intercept_d[[i]],
                                     legend = legend_text[[i]], legend_x = "topright",
                                     legend_col = 1, xaxt = xaxt[i], las = 1))
mtext(text = y_labels, side = 2, line = 3, at = 1.05 * ylim[2])
dev.off()

############# Plot year effect for c_age (Figure E.15 in appendix) #############

# plot on clr-level arranged as matrix vertical
le <- 8
re <- 1.5
top_bot <- 5
mar <- c(list(c(0, le, top_bot, re), c(0, le, 0, re), c(top_bot, le, 0, re)))
y_labels <- list("", expression(paste(clr, " [", hat(g)[c_age](year), "]")), "")
xaxt <- c("n", "n", "s")
year_cgroup <- list(effects$year_1, effects$year_2, effects$year_3)
legend_text <- paste0("c_age = ", c("0-6", "7-18", "other"))
ylim_clr_year <- range(year_cgroup, na.rm = TRUE)

pdf(paste0(path, "estimated_year_c_age_clr.pdf"), width = 5, height = 6.9)
layout(matrix(1:3, ncol = 1), heights = c(1, 0.74, 1))
sapply(1:3, function(i) plot_effects(year_cgroup[[i]], col = color, lty = lty_year, pdf = FALSE,
                                     ylim = ylim_clr_year, mar = mar[[i]],
                                     legend = legend_text[[i]], legend_x = "bottom",
                                     legend_col = 1, legend_cex = 1.5, # y_at = -3:3,
                                     ylab = "", cex.axis = 1.5, cex.lab = 1.5,
                                     xaxt = xaxt[i], las = 1, abline_h = NULL))
mtext(text = y_labels[[2]], side = 2, line = 4.5, at = max(ylim_clr_year) + 0.55 * diff(range(ylim_clr_year)))
dev.off()

# Plot Intercept \oplus year \oplus year_c_age on Bayes level
y_labels <- expression(paste("", "", hat(paste("", beta)), phantom()[{paste("0")}], "", paste(" "), "",
                             symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age")}], paste(" "), "",
                             symbol("\305"), paste(" "), "", hat(paste("", g)),  "", (year), paste(" "), "",
                             symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age")}], "", (year)))
year_cgroup_intercept <- list(effects$year_1_intercept_pdf,
                              effects$year_2_intercept_pdf,
                              effects$year_3_intercept_pdf)
year_cgroup_intercept_d <- list(effects$year_1_intercept_d_pdf,
                                effects$year_2_intercept_d_pdf,
                                effects$year_3_intercept_d_pdf)
ylim <- c(0, round(max(unlist(lapply(year_cgroup_intercept, max, na.rm = TRUE))), 1))

pdf(paste0(path, "estimated_year_c_age_intercept.pdf"), width = 5, height = 6.9)
layout(matrix(1:3, ncol = 1), heights = c(1, 0.74, 1))
sapply(1:3, function(i) plot_effects(year_cgroup_intercept[[i]], col = color, lty = lty_year, pdf = FALSE,
                                     ylim = ylim, mar = mar[[i]], # lwd = 1,
                                     effect_matrix_d = year_cgroup_intercept_d[[i]], pch = 1,
                                     legend = legend_text[[i]], legend_x = "topright",
                                     legend_col = 1, legend_cex = 1.5,
                                     ylab = "", cex.axis = 1.5, cex.lab = 1.5,
                                     xaxt = xaxt[i], las = 1))
mtext(text = y_labels, side = 2, line = 4.5, at = max(ylim) + max(ylim) / 2, cex = 1)
dev.off()

############## Plot year effects for interaction West/East and c_age #############
##############              (Figure E.16 in appendix)              #############

# plot on clr-level arranged as matrix vertical
le <- 8
re <- 1.5
top_bot <- 5
mar <- c(list(c(0, le, top_bot, re), c(0, le, 0, re), c(top_bot, le, 0, re)))

year_WestEast_cgroup <- list(effects$year_West_1, effects$year_West_2, effects$year_West_3,
                           effects$year_East_1, effects$year_East_2, effects$year_East_3)
year_WestEast_cgroup_d <- list(effects$year_West_1_d, effects$year_West_2_d, effects$year_West_3_d,
                             effects$year_East_1_d, effects$year_East_2_d, effects$year_East_3_d)
legend_text <- paste0("c_age = ", c("0-6", "7-18", "other"))
y_labels_West <- expression(paste(clr, " [", hat(g)["c_age, West"](year), "]"))
y_labels_East <- expression(paste(clr, " [", hat(g)["c_age, East"](year), "]"))
xaxt <- c(rep("n", 2), "s")
ylim_clr_year_West_East_c_age <- range(sapply(year_WestEast_cgroup, range, na.rm = TRUE))
at_clr_year_West_East_c_age <- max(ylim_clr_year_West_East_c_age) + 0.55 * diff(range(ylim_clr_year_West_East_c_age))

pdf(paste0(path, "estimated_year_West_c_age_clr.pdf"), width = 5, height = 6.7)
layout(matrix(1:3, ncol = 1), heights = c(1, 0.74, 1))
sapply(1:3, function(i) plot_effects(year_WestEast_cgroup[[i]], col = color,
                                     lty = lty_year, pdf = FALSE,
                                     mar = mar[[i]], ylim = ylim_clr_year_West_East_c_age,
                                     legend = legend_text[[i]], legend_x = "top",
                                     legend_col = 1, legend_cex = 1.5,
                                     # y_at = seq(-1.5, 0.5, by = 0.5),
                                     ylab = "", cex.axis = 1.5, cex.lab = 1.5,
                                     xaxt = xaxt[i], las = 1, abline_h = NULL))
mtext(text = y_labels_West, side = 2, line = 4.5, at = at_clr_year_West_East_c_age, cex = 1)
dev.off()

pdf(paste0(path, "estimated_year_East_c_age_clr.pdf"), width = 5, height = 6.7)
layout(matrix(1:3, ncol = 1), heights = c(1, 0.74, 1))
sapply(4:6, function(i) plot_effects(year_WestEast_cgroup[[i]], col = color,
                                     lty = lty_year, pdf = FALSE,
                                     mar = mar[[i-3]], ylim = ylim_clr_year_West_East_c_age,
                                     legend = legend_text[[i-3]], legend_x = "top",
                                     legend_col = 1, legend_cex = 1.5,
                                     # y_at = seq(-1.5, 0.5, by = 0.5),
                                     ylab = "", cex.axis = 1.5, cex.lab = 1.5,
                                     xaxt = xaxt[i-3], las = 1, abline_h = NULL))
mtext(text = y_labels_East, side = 2, line = 4.5, at = at_clr_year_West_East_c_age, cex = 1)
dev.off()

# Plot Intercept + West_East + c_age + West_East_c_age + year_intercept + year_West_East + year_cgroup + year_West_East_cgroup on Bayes level
y_labels_West <- expression(paste("", "", hat(paste("", beta)), phantom()[{paste("0")}], "", paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "West")}], paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age")}], paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age, West")}], paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)),  "", (year), paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "West")}], "", (year), paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age")}], "", (year), paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age, West")}], "", (year)))
y_labels_East <- expression(paste("", "", hat(paste("", beta)), phantom()[{paste("0")}], "", paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "East")}], paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age")}], paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", beta)), phantom()[{paste("", "c_age, East")}], paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)), "", (year), paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "East")}], "", (year), paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age")}], "", (year), paste(" "), "",
                                 symbol("\305"), paste(" "), "", hat(paste("", g)), phantom()[{paste("", "c_age, East")}], "", (year)))

year_West_cgroup_intercept <- list(effects$year_West_1_intercept_pdf, effects$year_West_2_intercept_pdf, effects$year_West_3_intercept_pdf)
year_East_cgroup_intercept <- list(effects$year_East_1_intercept_pdf, effects$year_East_2_intercept_pdf, effects$year_East_3_intercept_pdf)
year_West_cgroup_intercept_d <- list(effects$year_West_1_intercept_d_pdf, effects$year_West_2_intercept_d_pdf, effects$year_West_3_intercept_d_pdf)
year_East_cgroup_intercept_d <- list(effects$year_East_1_intercept_d_pdf, effects$year_East_2_intercept_d_pdf, effects$year_East_3_intercept_d_pdf)

ylim <- c(0, round(max(unlist(lapply(year_West_cgroup_intercept, max, na.rm = TRUE)),
                       unlist(lapply(year_East_cgroup_intercept, max, na.rm = TRUE))), 1))

pdf(paste0(path, "estimated_year_West_c_age_intercept.pdf"), width = 5, height = 6.7)
layout(matrix(1:3, ncol = 1), heights = c(1, 0.74, 1))
sapply(1:3, function(i) plot_effects(year_West_cgroup_intercept[[i]], col = color,
                                     lty = lty_year, pdf = FALSE,
                                     ylim = ylim, mar = mar[[i]], pch = 1,
                                     effect_matrix_d = year_West_cgroup_intercept_d[[i]],
                                     legend = legend_text[[i]], legend_x = "topright",
                                     legend_col = 1, legend_cex = 1.5,
                                     y_at = seq(0, 2.5, by = 0.5),
                                     ylab = "", cex.axis = 1.5, cex.lab = 1.5,
                                     xaxt = xaxt[i], las = 1))
mtext(text = y_labels_West, side = 2, line = 4.5, at = 1.57 * ylim[2], cex = 1)
dev.off()

pdf(paste0(path, "estimated_year_East_c_age_intercept.pdf"), width = 5, height = 6.7)
layout(matrix(1:3, ncol = 1), heights = c(1, 0.74, 1))
sapply(1:3, function(i) plot_effects(year_East_cgroup_intercept[[i]], col = color,
                                     lty = lty_year, pdf = FALSE,
                                     ylim = ylim, mar = mar[[i]], pch = 1,
                                     effect_matrix_d = year_East_cgroup_intercept_d[[i]],
                                     legend = legend_text[[i]], legend_x = "topright",
                                     legend_col = 1, legend_cex = 1.5,
                                     y_at = seq(0, 2.5, by = 0.5),
                                     ylab = "", cex.axis = 1.5, cex.lab = 1.5,
                                     xaxt = xaxt[i], las = 1))
mtext(text = y_labels_East, side = 2, line = 4.5, at = 1.57 * ylim[2], cex = 1)
dev.off()

par(def.par)
rm(list = ls())

