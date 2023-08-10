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
                      levels = c("old", "new"))
dta$region <- factor(dta$region)
dta$c_age <- factor(dta$child_group, levels = c(3, 2, 1), ordered = TRUE)
dta$child_group <- NULL

step_size <- 0.01
n_bins <- 1/ step_size

# source("./experiments/data_prep.R")

dta_grouped <- dta[, group_id := .GRP,
                   keyby = .(old_new # region,
                             #c_age,
                             #syear
                   )][, .(share,
                          old_new,
                          # region,
                          # c_age,
                          # syear,
                          group_id)]

dta_est <- construct_dataset(dta_grouped, response = list("share"),
                             step_size = step_size)

# weights for discrete points are 1, equal equidistant stepsize for all covariate
# combinations (histograms) in continuous component
dta_est$Delta <- ifelse(dta_est$share %in% 0:1, 1, step_size)

# Comparison ordered vs. non-ordered factor (ordered should be reference coding, non-ordered dummy-coding?)
model_old_new_dummy <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                                       mc = FALSE, np = FALSE) # default: no penalty for discrete component
                           + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                                mc = FALSE, np = FALSE, by = old_new)
                           + as.factor(group_id) - 1 + offset(log(Delta)), # no scalar intercept, but one intercept per covariate combination
                           data = dta_est, method = "REML", family = poisson())

dta$old_new <- ordered(dta$old_new)
dta_grouped <- dta[, group_id := .GRP,
                   keyby = .(old_new # region,
                             #c_age,
                             #syear
                   )][, .(share,
                          old_new,
                          # region,
                          # c_age,
                          # syear,
                          group_id)]

dta_est <- construct_dataset(dta_grouped, response = list("share"),
                             step_size = step_size)
dta_est$Delta <- ifelse(dta_est$share %in% 0:1, 1, step_size)

model_old_new <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                                 mc = FALSE, np = FALSE) # default: no penalty for discrete component
                     + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                          mc = FALSE, np = FALSE, by = old_new)
                     # + ti(share, syear, bs = c("md","ps"), m = list(c(2,2), c(2,2)),
                     #      k = c(12, 8), mc = c(FALSE, TRUE), np = FALSE) # k = c(54, 16)
                     + as.factor(group_id) - 1 + offset(log(Delta)), # no scalar intercept, but one intercept per covariate combination
                     data = dta_est, method = "REML", family = poisson())

######### Plot predictions on clr- and Bayes-level

X_old_new_dummy <- model.matrix(model_old_new_dummy)
# remove intercepts per covariate combination (not of interest for estimated densities)
intercepts <- which(grepl("group_id", colnames(X_old_new_dummy)))
X_old_new_dummy <- X_old_new_dummy[, -intercepts]
theta_hat_old_new_dummy <- model_old_new_dummy$coefficients[-intercepts]
f_hat_clr_old_new_dummy <- X_old_new_dummy %*% theta_hat_old_new_dummy

share <- unique(dta_est$share)
old_new <- unique(dta_est$old_new)
f_hat_clr_old_new_dummy <- matrix(f_hat_clr_old_new_dummy, nrow = length(share),
                            ncol = length(old_new))

X_old_new <- model.matrix(model_old_new)
# remove intercepts per covariate combination (not of interest for estimated densities)
intercepts <- which(grepl("group_id", colnames(X_old_new)))
X_old_new <- X_old_new[, -intercepts]
theta_hat_old_new <- model_old_new$coefficients[-intercepts]
f_hat_clr_old_new <- X_old_new %*% theta_hat_old_new

f_hat_clr_old_new <- matrix(f_hat_clr_old_new, nrow = length(share),
                                  ncol = length(old_new))


# Load (amongst others) function to plot densities
source("./experiments/help_functions.R")

library(RColorBrewer)
old_new_col <- brewer.pal(3, "Set1")[2:1]

### Plot clr-level:
par(mfrow = c(1, 2))
matplot(share[2:101], f_hat_clr_old_new_dummy[2:101,], type = "l", ylab = "estimated clr-density",
        xlab = "woman's income share", ylim = range(f_hat_clr_old_new_dummy), col = old_new_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_clr_old_new_dummy[c(1, 102), ]), w = 0.02, col = old_new_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

matplot(share[2:101], f_hat_clr_old_new[2:101,], type = "l", ylab = "estimated clr-density",
        xlab = "woman's income share", ylim = range(f_hat_clr_old_new), col = old_new_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_clr_old_new[c(1, 102), ]), w = 0.02, col = old_new_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

mean(abs(f_hat_clr_old_new_dummy - f_hat_clr_old_new))
max(abs(f_hat_clr_old_new_dummy - f_hat_clr_old_new))

# Transform to Bayes-level
f_hat_old_new_dummy <- clr_inv(f_hat_clr_old_new_dummy)
f_hat_old_new <- clr_inv(f_hat_clr_old_new)

### Plot Bayes-level:
matplot(share[2:101], f_hat_old_new_dummy[2:101,], type = "l", ylab = "estimated density",
        xlab = "woman's income share", ylim = c(0, max(f_hat_old_new_dummy)),
        col = old_new_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_old_new_dummy[c(1, 102), ]), w = 0.02, col = old_new_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

matplot(share[2:101], f_hat_old_new[2:101,], type = "l", ylab = "estimated density",
        xlab = "woman's income share", ylim = c(0, max(f_hat_old_new)), col = old_new_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_old_new[c(1, 102), ]), w = 0.02, col = old_new_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

mean(abs(f_hat_old_new_dummy - f_hat_old_new))
max(abs(f_hat_old_new_dummy - f_hat_old_new))

# Not exactly the same, but pretty close (and difference not visible)

################### Include interaction old_new/c_age ####################

# ### Is it possible to specify more than one by-variable?
# dta_grouped <- dta[, group_id := .GRP,
#                    keyby = .(old_new, # region,
#                              c_age
#                              #syear
#                    )][, .(share,
#                           old_new,
#                           # region,
#                           c_age,
#                           # syear,
#                           group_id)]
#
# dta_est <- construct_dataset(dta_grouped, response = list("share"),
#                              step_size = step_size)
# dta_est$Delta <- ifelse(dta_est$share %in% 0:1, 1, step_size)
#
#
# model_interaction <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
#                                      mc = FALSE, np = FALSE) # default: no penalty for discrete component
#                          + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
#                               mc = FALSE, np = FALSE, by = old_new)
#                          + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
#                               mc = FALSE, np = FALSE, by = c_age)
#                          + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
#                               mc = FALSE, np = FALSE, by = c(old_new, c_age))
#                          # + ti(share, syear, bs = c("md","ps"), m = list(c(2,2), c(2,2)),
#                          #      k = c(12, 8), mc = c(FALSE, TRUE), np = FALSE) # k = c(54, 16)
#                          + as.factor(group_id) - 1 + offset(log(Delta)), # no scalar intercept, but one intercept per covariate combination
#                          data = dta_est, method = "REML", family = poisson())
#
# # No... I also tried by = list(old_new, c_age), which also gave an error...

### Create new factor variable for interaction (like in 1st paper)
dta$old_new_c_age <- paste0(dta$old_new, "_", dta$c_age)
dta$old_new_c_age <- factor(dta$old_new_c_age, levels = c(paste0("old_", 3:1), paste0("new_", 3:1)),
                                  ordered = TRUE)

dta_grouped <- dta[, group_id := .GRP,
                   keyby = .(old_new, # region,
                             c_age,
                             old_new_c_age
                             #syear
                   )][, .(share,
                          old_new,
                          # region,
                          c_age,
                          old_new_c_age,
                          # syear,
                          group_id)]

dta_est <- construct_dataset(dta_grouped, response = list("share"),
                             step_size = step_size)
dta_est$Delta <- ifelse(dta_est$share %in% 0:1, 1, step_size)

model_interaction <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                                      mc = FALSE, np = FALSE) # default: no penalty for discrete component
                          + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                               mc = FALSE, np = FALSE, by = old_new)
                          + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                               mc = FALSE, np = FALSE, by = c_age)
                          + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                               mc = FALSE, np = FALSE, by = old_new_c_age)
                          + as.factor(group_id) - 1 + offset(log(Delta)), # no scalar intercept, but one intercept per covariate combination
                          data = dta_est, method = "REML", family = poisson())

######### Plot predictions on clr- and Bayes-level

X_interaction <- model.matrix(model_interaction)
# remove intercepts per covariate combination (not of interest for estimated densities)
intercepts <- which(grepl("group_id", colnames(X_interaction)))
X_interaction <- X_interaction[, -intercepts]
theta_hat_interaction <- model_interaction$coefficients[-intercepts]
f_hat_clr_interaction <- X_interaction %*% theta_hat_interaction

share <- unique(dta_est$share)
old_new_c_age <- unique(dta_est$old_new_c_age)
f_hat_clr_interaction <- matrix(f_hat_clr_interaction, nrow = length(share),
                                  ncol = length(old_new_c_age))

# Load (amongst others) function to plot densities
source("./experiments/help_functions.R")

# library(RColorBrewer)
# old_new_col <- brewer.pal(3, "Set1")[2:1]
c_age_col <- brewer.pal(n = 3, "Dark2")[3:1]
interaction_col <- c(brewer.pal(n = 3, "Dark2")[3:1], brewer.pal(n = 3, "Accent")[c(2, 3, 1)])

### Plot clr-level:
matplot(share[2:101], f_hat_clr_interaction[2:101,], type = "l", ylab = "estimated clr-density",
        xlab = "woman's income share", ylim = range(f_hat_clr_interaction), col = interaction_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_clr_interaction[c(1, 102), ]), w = 0.02, col = interaction_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)

# Transform to Bayes-level
f_hat_interaction <- clr_inv(f_hat_clr_interaction)

### Plot Bayes-level:
matplot(share[2:101], f_hat_interaction[2:101,], type = "l", ylab = "estimated density",
        xlab = "woman's income share", ylim = c(0, max(f_hat_interaction)),
        col = interaction_col, lty = 1)
my_minus(c(-0.005, 1.005), t(f_hat_interaction[c(1, 102), ]), w = 0.02, col = interaction_col, lwd = 0.9)
abline(h = 0, col = "grey", lty = 2)


# predict(..., type = "terms") contains estimated effects on clr-level with columns:
# - 1: intercepts per covariate combination (not of interest for us)
# - 2: smooth intercept (reference old, no children)
# - 3: smooth effect for new
# - 4-5: smooth effects for child groups 2 and 1
# - 6-10: smooth effects for interaction of old_new and c_age
pred_terms <- predict(model_interaction, type = "terms")

intercept <- pred_terms[1:102, 2]


### Old_new effect (old: reference)
estimated_old_new_clr <- matrix(pred_terms[, 3], nrow = length(share),
                                ncol = length(old_new_c_age))[, c(1, 4)]
# column 1-3 (old_3, old_2, old_1) and 4-6 (new_3, new_2, new_1) are identical,
# respectively, with each of the colums  1-3 containing the old-effect and each
# of the columns 4-6 containing the new-effect

# Plot on clr-level
plot_effects(t(estimated_old_new_clr), col = old_new_col, pdf = FALSE, # path = path,
             abline_h = NULL, ylab = expression(paste("clr [",hat(beta)[West_East], "]")),
             legend = c("West_East = West", "West_East = East"), las = 1,
             legend_x = "bottomright", legend_cex = 0.7)

# Intercept + old_new
estimated_old_new_intercept_clr <- apply(estimated_old_new_clr, 2, "+", intercept)

# Transform to Bayes-level
estimated_old_new_intercept <- clr_inv(estimated_old_new_intercept_clr)
estimated_old_new_intercept_d <- get_discrete_pdf(estimated_old_new_intercept, per_col = TRUE)

# symbol("\305") corresponds to \oplus
ylab_intercept_old_new <- expression(paste("", "", hat(paste("", beta)),
                                           phantom()[{paste("0")}], "", paste(" "), "",
                                           symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                           phantom()[{paste("", "West_East")}], ""))


plot_effects(t(estimated_old_new_intercept), col = old_new_col, pdf = FALSE, #path = path,
             pch = 1, effect_matrix_d = t(estimated_old_new_intercept_d),
             ylab = ylab_intercept_old_new, las = 1, legend_cex = 0.7,
             legend = c("West_East = West", "West_East = East"), legend_x = "topright")


### c_age effect (old: reference)
estimated_c_age_clr <- Reduce("+", lapply(4:5, function(k) matrix(pred_terms[, k], nrow = length(share),
                                                                        ncol = length(old_new_c_age))[, 1:3]))

# Intercept + c_age
estimated_c_age_intercept_clr <- apply(estimated_c_age_clr, 2, "+", intercept)
# Transform to Bayes-level
estimated_c_age_intercept <- clr_inv(estimated_c_age_intercept_clr)

# c_age_col <- brewer.pal(n = 3, "Dark2")[3:1]
c_age_legend <- c("c_age = other", "c_age = 0-6", "c_age = 7-18")

# Plot c_age on clr-level
plot_effects(t(estimated_c_age_clr), col = c_age_col, pdf = FALSE, #path = path,
             abline_h = NULL, ylab = expression(paste(clr, " [", hat(beta)["c_age"], "]")),
             las = 1, legend = c_age_legend, legend_x = "bottomleft", legend_cex = 0.7)

# Plot intercept \oplus c_age on Bayes level
ylab_intercept_cgroup <- expression(paste("", "", hat(paste("", beta)),
                                          phantom()[{paste("0")}], "", paste(" "), "",
                                          symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                          phantom()[{paste("", "c_age")}], ""))

plot_effects(t(estimated_c_age_intercept), col = c_age_col, pdf = FALSE, #path = path,
             # effect_matrix_d = effects$c_age_intercept_d_pdf, pch = 1,
             ylab = ylab_intercept_cgroup, las = 1, legend_cex = 0.7,
             legend = c_age_legend, legend_x = "topright")

### interaction effect old_new and c_age (old, 3: reference)
estimated_old_new_c_age_clr <- Reduce("+", lapply(6:10, function(k) matrix(pred_terms[, k], nrow = length(share),
                                                                  ncol = length(old_new_c_age))))

interaction_legend <- paste0("c_age = ", c("0-6", "7-18", "other"),
                             ", West_East = ", c(rep("West", 3), rep("East", 3)))

# Plot interaction on clr-level
plot_effects(t(estimated_old_new_c_age_clr), col = interaction_col, pdf = FALSE, # path = path,
             ylab = expression(paste(clr, " [", hat(beta)["c_age, West_East"], "]")),
             legend = interaction_legend, legend_x = "bottom", legend_ncol = 2,
             legend_col = interaction_col[c(3:1, 6:4)], legend_cex = 0.7, abline_h = NULL, las = 1)

# Intercept + old_new_c_age
estimated_old_new_c_age_intercept_clr <- apply(estimated_old_new_c_age_clr, 2, "+", intercept)
# Transform to Bayes-level
estimated_old_new_c_age_intercept <- clr_inv(estimated_old_new_c_age_intercept_clr)

# Interactions with East (new) seam a little bit weird for non-reference groups:
# East-c_age2 (7-18) is straight line, East-c_age1 (0-6) is constantly (almost)
# zero... (both clr-level)
max(abs(pred_terms[, 10]))

# Plot intercept \oplus old_new \oplus c_age \oplus old_new_c_age on Bayes level
ylab_old_new_c_age_1 <- expression(paste("", "", hat(paste("", beta)),
                                         phantom()[{paste("0")}], "", paste(" "), "",
                                         symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                         phantom()[{paste("", "West_East")}], "", paste(" "), "",
                                         symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                         phantom()[{paste("", "c_age")}], ""))
ylab_old_new_c_age_2 <- expression(paste(symbol("\305"), paste(" "), "", hat(paste("", beta)),
                                         phantom()[{paste("", "c_age, West_East")}], ""))
at_old_new_c_age <- max(estimated_old_new_c_age_intercept, na.rm = TRUE)/2


plot_effects(t(estimated_old_new_c_age_intercept), col = interaction_col, # path = path,
             pdf = FALSE,  pch = 1, ylab = "",
             las = 1, # effect_matrix_d = effects$old_new_c_age_intercept_d_pdf,
             legend_col = interaction_col[c(3:1, 6:4)],
             legend_cex = 0.7, legend = interaction_legend, legend_x = "topright")
mtext(text = ylab_old_new_c_age_1, side = 2, line = 3.8, at = at_old_new_c_age)
mtext(text = ylab_old_new_c_age_2, side = 2, line = 2.5, at = at_old_new_c_age)

