################################ Initialization ################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/.."))

library(data.table)
library(mgcv)
source("./R/mixed_density_smooth.R")

# read in data to data.table
dta <- as.data.table(readRDS("../income_share_data.rds"))

step_size <- 0.01
n_bins <- 1/ step_size

########################## Univariate (no covariates) ##########################

# create bins and get bin mid points and counts
grid_hist <- seq(from = 0, to = 1, by = step_size)
hist_dens <- hist(dta[share > 0 & share < 1]$share, breaks = grid_hist,
                  right = FALSE) # In the paper, we defined the binning intervals to be left-closed
data_d <-dta[!(share > 0 & share < 1)][, .N, by=share]
colnames(data_d) <- c("share", "counts")

counts_dens <- hist_dens$counts
mids_dens <- hist_dens$mids
dta_dens_c <- as.data.table(cbind(counts_dens, mids_dens))
colnames(dta_dens_c) <- c("counts", "share")
dta_dens <- rbind(dta_dens_c, data_d)[order(share)]

# bin widths (continuous component) and integration weights (discrete component)
# required as (log-)offset in model specification
Delta <- c(1, rep(step_size, length(mids_dens)), 1)

# Notes on model specification:
# - bs = "md" includes the integrate-to-zero constraint already and no additional
# centering shall be applied -> use ti(..., mc = FALSE), since s() and te() don't
# have an option to turn centering constraints off!
# - np = TRUE (default) implements a reparametrization, which is not desired in
# our case -> use ti(..., np = FALSE)
# - offset(log(Delta)) adds the necessary additive term in the predictor that
# includes binwidths/dirac weights into the estimation
model_null <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                              mc = FALSE, np = FALSE) + # default: no penalty for discrete component
                    offset(log(Delta)),
                  data = dta_dens, method = "REML", family = poisson())

model_ridge <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                               xt = list(penalty_discrete = 0), # include penalty for discrete component
                               mc = FALSE, np = FALSE) + offset(log(Delta)),
                   data = dta_dens, method = "REML", family = poisson())
# model_null$coefficients == model_ridge$coefficients
# model_null$smooth[[1]]$margin[[1]]$knots

# estimated clr-densities (exclude intercepts of covariate combinations, here
# only one):
X_null <- model.matrix(model_null)[, -1]
theta_hat_null <- model_null$coefficients[-1]
f_hat_clr_null <- c(X_null %*% theta_hat_null)

X_ridge <- model.matrix(model_ridge)[, -1]
theta_hat_ridge <- model_ridge$coefficients[-1]
f_hat_clr_ridge <- c(X_ridge %*% theta_hat_ridge)

# Plot on clr-level:
plot(mids_dens, f_hat_clr_null[2:101], type = "l", ylab = "estimated clr-density",
     xlab = "woman's income share", ylim = range(f_hat_clr_null, f_hat_clr_ridge))
points(c(0, 1), f_hat_clr_null[c(1, 102)])
abline(h = 0, col = "grey", lty = 2)
lines(mids_dens, f_hat_clr_ridge[2:101], col = 2, lty = 2)
points(c(0, 1), f_hat_clr_ridge[c(1, 102)], col = 2, pch = 4)
legend("bottomleft", legend = c("No penalty", "Ridge penalty"), lty = 1:2,
       col = 1:2, bty = "n", title = "Penalty for discrete component:")

# Since both estimates are very similar, we transform only one model to Bayes-level
f_hat_null <- exp(f_hat_clr_null) / sum(exp(f_hat_clr_null[c(1, 102)]), mean(exp(f_hat_clr_null[2:101])))

# Plot histogram on relative (i.e., density) scale and add estimated density
zeros_dens <- length(which(dta$share == 0)) / nrow(dta)
dual_earner_dens <- sum(counts_dens) / nrow(dta) * counts_dens / (step_size * sum(counts_dens))
ones_dens <- length(which(dta$share == 1)) / nrow(dta)

barplot(dual_earner_dens, width = diff(mids_dens)[1], space = 0,
        main = "Histogram of SOEP data with estimated density",
        ylim = c(0, max(zeros_dens, dual_earner_dens, ones_dens, f_hat_null) + 0.01))
axis(1)
points(c(0, 1), c(zeros_dens, ones_dens))
lines(mids_dens, f_hat_null[2:101], type = "l", ylim = range(f_hat_null), col = "red")
points(c(0, 1), f_hat_null[c(1, 102)], col = "red", pch = 4)


# Plots on Response-level of Poisson model
library(ggplot2)

dens_unnorm <- predict(model_null, type = "response")
gg_data_uni <- data.table(share = dta_dens$share, predictions = dens_unnorm)

# plot counts underlying the Poisson model ("histogram") and add estimated function
ggplot(dta, aes(x = share)) +
  geom_histogram(bins = n_bins) +
  labs(x = "share of woman's income", y = "number of observations") +
  geom_line(data = gg_data_uni,
            aes(x = share, y = predictions),
            color = "red",
            linewidth = 1)

