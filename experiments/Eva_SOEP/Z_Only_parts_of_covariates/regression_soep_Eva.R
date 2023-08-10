current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/.."))

## This code implements the poisson density regression approach for the SOEP dataset
library(data.table)
library(mgcv)

# source("./experiments/data_prep.R")
# source("./experiments/mixed_density_smooth.R")
source("./R/mixed_density_smooth.R")
source("./experiments/data_prep.R")

# read in data to data.table
dta <- as.data.table(readRDS("../income_share_data.rds"))
dta$old_new <- ifelse(dta$region %in% c("northeast", "east"), "new", "old")
dta$old_new <- factor(ifelse(dta$region %in% c("northeast", "east"), "new", "old"),
                      levels = c("old", "new"), ordered = TRUE) # ordered factor for reference coding later
dta$region <- factor(dta$region)
dta$child_group <- factor(dta$child_group)

step_size <- 0.01
n_bins <- 1/ step_size

############################### Model with year ################################
dta_syear <- dta[, group_id := .GRP,
                   keyby = .(#region,
                     #child_group,
                     syear)][, .(share,
                                 # old_new,
                                 # region,
                                 # child_group,
                                 syear,
                                 group_id)]

dta_syear <- construct_dataset(dta_syear, response = list("share"),
                             step_size = step_size)

# weights for discrete points are 1, equal equidistant stepsize for all covariate
# combinations (histograms) in continuous component
dta_syear$Delta <- ifelse(dta_syear$share %in% 0:1, 1, step_size)

# Notes on model specification:
# - bs = "md" includes the integrate-to-zero constraint already and no additional
# centering shall be applied -> use ti(..., mc = FALSE), since s() and te() don't
# have an option to turn centering constraints off!
# - np = TRUE (default) implements a reparametrization, which is not desired in
# our case -> use ti(..., np = FALSE)
# - offset(log(Delta)) adds the necessary additive term in the predictor that
# includes binwidths/dirac weights into the estimation
model_syear <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12, # k = 54 would equal number of knots in first paper, however, there there were smooth densities estimated first; here, no pre-smoothing happens and so many knots lead to overfitting!
                               mc = FALSE, np = FALSE) # default: no penalty for discrete component
                   + ti(share, syear, bs = c("md","ps"), m = list(c(2,2), c(2,2)),
                        k = c(12, 8), mc = c(FALSE, TRUE), np = FALSE) # k = c(54, 16)
                   + as.factor(group_id) - 1 + offset(log(Delta)), # no scalar intercept, but one intercept per covariate combination
                   data = dta_syear, method = "REML", family = poisson())


######### Plot predictions on clr- and Bayes-level

X_syear <- model.matrix(model_syear)
# remove intercepts per covariate combination (not of interest for estimated densities)
intercepts_syear <- which(grepl("group_id", colnames(X_syear)))
X_syear <- X_syear[, -intercepts_syear]
theta_hat_syear <- model_syear$coefficients[-intercepts_syear]
f_hat_clr_syear <- X_syear %*% theta_hat_syear

share <- unique(dta_est$share)
year <- unique(dta_est$syear)
f_hat_clr_syear <- matrix(f_hat_clr_syear, nrow = length(share),
                          ncol = length(year))

# Load (amongst others) function to plot densities
source("./experiments/help_functions.R")

# Color gradient over the 33 observed years: from dark purple (1980s) to light
# purple (1990s), red/orange (2000s) to yellow (2016)
color <- viridis::plasma(33)
# Different line types for further differentiation
lty <- rep(c(1, 2, 4, 5), 9)

### Plot clr-level:
matplot(share[2:101], f_hat_clr_syear[2:101,], type = "l", ylab = "estimated clr-density",
        xlab = "woman's income share", ylim = range(f_hat_clr_syear), col = color, lty = lty)
my_minus(c(-0.005, 1.005), t(f_hat_clr_syear[c(1, 102), ]), w = 0.02, col = color, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

# Transform to Bayes-level
f_hat_syear <- clr_inv(f_hat_clr_syear)

### Plot Bayes-level:
matplot(share[2:101], f_hat_syear[2:101,], type = "l", ylab = "estimated density",
        xlab = "woman's income share", ylim = c(0, max(f_hat_syear)), col = color, lty = lty)
my_minus(c(-0.005, 1.005), t(f_hat_syear[c(1, 102), ]), w = 0.02, col = color, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

# For each year: Plot histogram on relative (i.e., density) scale and add estimated density
zeros_dens <- lapply(unique(dta_est$group_id),
                     function(g) dta_est[group_id == g & share == 0]$counts / nrow(dta_grouped[group_id == g]))
dual_earner_dens <- lapply(unique(dta_est$group_id),
                           function(g) {
                             counts_dens <- dta_est[group_id == g & !(share %in% 0:1)]$counts
                             sum(counts_dens) / nrow(dta_grouped[group_id == g]) * counts_dens / (step_size * sum(counts_dens))
                           })
ones_dens <- lapply(unique(dta_est$group_id),
                    function(g) dta_est[group_id == g & share == 1]$counts / nrow(dta_grouped[group_id == g]))
# # check
# sapply(seq_len(length(year)), function(i) sum(c(zeros_dens[[i]],
#                                                 dual_earner_dens[[i]] * step_size,
#                                                 ones_dens[[i]])))
# apply(f_hat_syear, 2, integrate_num)
# apply(f_hat_clr_syear, 2, integrate_num)

library(manipulate)
manipulate({
  barplot(dual_earner_dens[[i]], width = max(diff(share)), space = 0,
          main = paste0("Year: ", year[i]),
          ylim = c(0, max(zeros_dens[[i]], dual_earner_dens[[i]], ones_dens[[i]], f_hat_syear[, i]) + 0.01))
  axis(1)
  points(c(0, 1), c(zeros_dens[[i]], ones_dens[[i]]))
  lines(share[2:101], f_hat_syear[2:101, i], type = "l", col = color[i])
  points(c(0, 1), f_hat_syear[c(1, 102), i], col = color[i], pch = 4)
}, i = slider(1, length(year)))

######### Plot estimated effects on clr- and Bayes-level

# predict(..., type = "terms") contains estimated effects on clr-level:
# - 1st column: intercepts per covariate combination (not of interest for us)
# - 2nd column: smooth intercept
# - 3rd column: smooth year effect
pred_terms_syear <- predict(model_syear, type = "terms")

# # Check:
# # 1st column cooresponds to intercepts per covariate combination
# sum(unique(pred_terms[, 1]) != model_syear$coefficients[1:length(year)])
# # sum of relevant effects (i.e., without covariate specific intercepts)
# # should be equal to predicted clr-densities
# test <- matrix(apply(pred_terms[, 2:3], 1, sum), nrow = length(share), ncol = length(year))
# sum(test != f_hat_clr_syear)
# max(abs(test - f_hat_clr_syear))

### smooth intercept

# # Old way to extract smooth intercept (without predict)
# smooth_int <- model_syear$smooth[[1]]
# X_syear_smooth_int <- smooth_int$margin[[1]]$X
# theta_hat_syear_smooth_int <- model_syear$coefficients[smooth_int$first.para:smooth_int$last.para]
# estimated_smooth_int_clr <- X_syear_smooth_int %*% theta_hat_syear_smooth_int
# # same 102 (intercept) values are repeated 33 times:
# # sum(estimated_smooth_int_clr[1:102] != estimated_smooth_int_clr[103:204])
# estimated_smooth_int_clr <- estimated_smooth_int_clr[1:102]
# sum(estimated_smooth_int_clr != pred_terms[1:102, 2])

# same 102 (intercept) values are repeated 33 times: sum(pred_terms[1:102, 2] != pred_terms[103:204, 2])
syear_estimated_smooth_int_clr <- unique(pred_terms_syear[, 2])
# Plot on clr-level
plot(share[2:101], syear_estimated_smooth_int_clr[2:101], type = "l", ylab = "estimated clr-intercept",
     xlab = "woman's income share")
my_minus(c(-0.005, 1.005), t(syear_estimated_smooth_int_clr[c(1, 102)]), w = 0.02, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

# Transform to Bayes-level
syear_estimated_smooth_int <- clr_inv(syear_estimated_smooth_int_clr)

### Plot Bayes-level:
plot(share[2:101], syear_estimated_smooth_int[2:101], type = "l", ylab = "estimated intercept",
     xlab = "woman's income share")
my_minus(c(-0.005, 1.005), t(syear_estimated_smooth_int[c(1, 102)]), w = 0.02, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

### Year effect
estimated_year_clr <- matrix(pred_terms[, 3], nrow = length(share),  ncol = length(year))

# Plot on clr-level
matplot(share[2:101], estimated_year_clr[2:101,], type = "l", ylab = "estimated clr-effect year",
        xlab = "woman's income share", col = color, lty = lty)
my_minus(c(-0.005, 1.005), t(estimated_year_clr[c(1, 102),]), w = 0.02, lwd = 0.9, col = color)
abline(h = 0, col = "grey", lty = 2)

# Transform to Bayes-level
estimated_year <- clr_inv(estimated_year_clr)

### Plot Bayes-level:
matplot(share[2:101], estimated_year[2:101,], type = "l", ylab = "estimated clr-effect year",
        xlab = "woman's income share", col = color, lty = lty)
my_minus(c(-0.005, 1.005), t(estimated_year[c(1, 102), ]), w = 0.02, col = color, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

### Check whether intercepts could be modeled smoothly in the covariate basis
# get counts per covariate compbinaiton, here: year
counts_per_year <- sapply(unique(dta_est$group_id),
                          function(g) nrow(dta_grouped[group_id == g]))
par(mfrow = c(1, 2))
plot(year, counts_per_year, ylab = expression(n^{(year)}),
     main = "Observaion numbers per year")
plot(year, model_syear$coefficients[1:length(year)],
     ylab = expression(hat(tau)^{(year)}), main = "Estimated intercepts")
# -> Does not look smooth, there are rather waves of count levels (compare Info-
# material 2_Intro_SOEP_de.pdf, slide 8 and 32) -> Maybe estimate one intercept per wave?


############################## Model with old_new ##############################

# source("./experiments/data_prep.R")

dta_grouped <- dta[, group_id := .GRP,
                   keyby = .(old_new # region,
                             #child_group,
                             #syear
                             )][, .(share,
                                         old_new,
                                         # region,
                                         # child_group,
                                         # syear,
                                         group_id)]

dta_est <- construct_dataset(dta_grouped, response = list("share"),
                             step_size = step_size)

# weights for discrete points are 1, equal equidistant stepsize for all covariate
# combinations (histograms) in continuous component
dta_est$Delta <- ifelse(dta_est$share %in% 0:1, 1, step_size)

# Notes on model specification:
# - bs = "md" includes the integrate-to-zero constraint already and no additional
# centering shall be applied -> use ti(..., mc = FALSE), since s() and te() don't
# have an option to turn centering constraints off!
# - np = TRUE (default) implements a reparametrization, which is not desired in
# our case -> use ti(..., np = FALSE)
# - offset(log(Delta)) adds the necessary additive term in the predictor that
# includes binwidths/dirac weights into the estimation
model_old_new <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                                 mc = FALSE, np = FALSE) # default: no penalty for discrete component
                     + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                          mc = FALSE, np = FALSE, by = old_new)
                     # + ti(share, syear, bs = c("md","ps"), m = list(c(2,2), c(2,2)),
                     #      k = c(12, 8), mc = c(FALSE, TRUE), np = FALSE) # k = c(54, 16)
                     + as.factor(group_id) - 1 + offset(log(Delta)), # no scalar intercept, but one intercept per covariate combination
                     data = dta_est, method = "REML", family = poisson())




######### Plot predictions on clr- and Bayes-level

X_old_new <- model.matrix(model_old_new)
# Design matrix columns corresponding to share are equal:
sum(X_old_new[1:102, 3:15] != X_old_new[103:204, 3:15])
# remove intercepts per covariate combination (not of interest for estimated densities)
intercepts <- which(grepl("group_id", colnames(X_old_new)))
X_old_new <- X_old_new[, -intercepts]
theta_hat_old_new <- model_old_new$coefficients[-intercepts]
f_hat_clr_old_new <- X_old_new %*% theta_hat_old_new

share <- unique(dta_est$share)
old_new <- unique(dta_est$old_new)
f_hat_clr_old_new <- matrix(f_hat_clr_old_new, nrow = length(share),
                          ncol = length(old_new))

# Load (amongst others) function to plot densities
source("./experiments/help_functions.R")

library(RColorBrewer)
old_new_col <- brewer.pal(3, "Set1")[2:1]

### Plot clr-level:
matplot(share[2:101], f_hat_clr_old_new[2:101,], type = "l", ylab = "estimated clr-density",
        xlab = "woman's income share", ylim = range(f_hat_clr_old_new), col = old_new_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_clr_old_new[c(1, 102), ]), w = 0.02, col = old_new_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

# Transform to Bayes-level
f_hat_old_new <- clr_inv(f_hat_clr_old_new)

### Plot Bayes-level:
matplot(share[2:101], f_hat_old_new[2:101,], type = "l", ylab = "estimated density",
        xlab = "woman's income share", ylim = c(0, max(f_hat_old_new)), col = old_new_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_old_new[c(1, 102), ]), w = 0.02, col = old_new_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)


# For each year: Plot histogram on relative (i.e., density) scale and add estimated density
zeros_dens <- lapply(unique(dta_est$group_id),
                     function(g) dta_est[group_id == g & share == 0]$counts / nrow(dta_grouped[group_id == g]))
dual_earner_dens <- lapply(unique(dta_est$group_id),
                           function(g) {
                             counts_dens <- dta_est[group_id == g & !(share %in% 0:1)]$counts
                             sum(counts_dens) / nrow(dta_grouped[group_id == g]) * counts_dens / (step_size * sum(counts_dens))
                           })
ones_dens <- lapply(unique(dta_est$group_id),
                    function(g) dta_est[group_id == g & share == 1]$counts / nrow(dta_grouped[group_id == g]))
# # check
# sapply(seq_len(length(old_new)), function(i) sum(c(zeros_dens[[i]],
#                                                 dual_earner_dens[[i]] * step_size,
#                                                 ones_dens[[i]])))
# apply(f_hat_old_new, 2, integrate_num)
# apply(f_hat_clr_old_new, 2, integrate_num)

library(manipulate)
manipulate({
  i <- ifelse(old_new == "old", 1, 2)
  barplot(dual_earner_dens[[i]], width = max(diff(share)), space = 0,
          main = old_new,
          ylim = c(0, max(zeros_dens[[i]], dual_earner_dens[[i]], ones_dens[[i]], f_hat_old_new[, i]) + 0.1))
  axis(1)
  points(c(0, 1), c(zeros_dens[[i]], ones_dens[[i]]))
  lines(share[2:101], f_hat_old_new[2:101, i], type = "l", col = old_new_col[i])
  points(c(0, 1), f_hat_old_new[c(1, 102), i], col = old_new_col[i], pch = 4)
  },
  old_new = picker("old", "new"))

######### Plot estimated effects on clr- and Bayes-level

# predict(..., type = "terms") contains estimated effects on clr-level:
# - 1st column: intercepts per covariate combination (not of interest for us)
# - 2nd column: smooth intercept
# - 3rd column: smooth old_new effect
pred_terms <- predict(model_old_new, type = "terms")

### smooth intercept

# same 102 (intercept) values are repeated 2 times: sum(pred_terms[1:102, 2] != pred_terms[103:204, 2])
estimated_smooth_int_clr <- unique(pred_terms[, 2])
# Plot on clr-level
plot(share[2:101], estimated_smooth_int_clr[2:101], type = "l", ylab = "estimated clr-intercept",
     xlab = "woman's income share")
my_minus(c(-0.005, 1.005), t(estimated_smooth_int_clr[c(1, 102)]), w = 0.02, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

# Transform to Bayes-level
estimated_smooth_int <- clr_inv(estimated_smooth_int_clr)

### Plot Bayes-level:
plot(share[2:101], estimated_smooth_int[2:101], type = "l", ylab = "estimated intercept",
     xlab = "woman's income share")
my_minus(c(-0.005, 1.005), t(estimated_smooth_int[c(1, 102)]), w = 0.02, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

### old_new effect
estimated_old_new_clr <- matrix(pred_terms[, 3], nrow = length(share),  ncol = length(old_new))

# Plot on clr-level
matplot(share[2:101], estimated_old_new_clr[2:101,], type = "l", ylab = "estimated clr-effect old_new",
        xlab = "woman's income share", col = old_new_col, lty = lty)
my_minus(c(-0.005, 1.005), t(estimated_old_new_clr[c(1, 102),]), w = 0.02, lwd = 0.9, col = old_new_col)
abline(h = 0, col = "grey", lty = 2)

# Transform to Bayes-level
estimated_old_new <- clr_inv(estimated_old_new_clr)

### Plot Bayes-level:
matplot(share[2:101], estimated_old_new[2:101,], type = "l", ylab = "estimated clr-effect old_new",
        xlab = "woman's income share", col = old_new_col, lty = lty)
my_minus(c(-0.005, 1.005), t(estimated_old_new[c(1, 102), ]), w = 0.02, col = old_new_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

# ### Check whether intercepts could be modeled smoothly in the covariate basis
# # get counts per covariate compbinaiton, here: old_new
# counts_per_old_new <- sapply(unique(dta_est$group_id),
#                           function(g) nrow(dta_grouped[group_id == g]))
# par(mfrow = c(1, 2))
# plot(old_new, counts_per_old_new, ylab = expression(n^{(old_new)}),
#      main = "Observaion numbers per old_new")
# plot(old_new, model_old_new$coefficients[1:length(old_new)],
#      ylab = expression(hat(tau)^{(old_new)}), main = "Estimated intercepts")

