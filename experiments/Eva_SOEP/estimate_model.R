################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/../.."))

# Package to deal with big datasets
library(data.table)
# Package to estimate GAMs (like our Poisson model)
library(mgcv)

# load self-written smoother for mixed reference measure
source("./R/mixed_density_smooth.R")
# load functions used to prepare data appropriately to be used in Poisson model
source("./experiments/data_prep.R")

# load data, add missing variables and transform characters to ordered factors
# (necessary for reference coding in gam(); lowest factor level is used as
# reference category)
dta <- as.data.table(readRDS("../income_share_data.rds"))
dta$West_East <- factor(ifelse(dta$region %in% c("northeast", "east"), "East", "West"),
                        levels = c("West", "East"), ordered = TRUE)
# dta$region <- factor(dta$region) # I don't know how to include region such that it is centered around West/East -> Need to talk to Almond
dta$c_age <- factor(dta$child_group, levels = c(3, 2, 1), ordered = TRUE)

# Create new factor variable for interaction of West_East and c_age
dta$West_East_c_age <- paste0(dta$West_East, "_", dta$c_age)
dta$West_East_c_age <- factor(dta$West_East_c_age, levels = c(paste0("West_", 3:1), paste0("East_", 3:1)),
                              ordered = TRUE)

# Normalize SOEP weight (compare to ?gam -> weights);
dta$weighting_factor_orig <- dta$weighting_factor
dta$weighting_factor <- dta$weighting_factor_orig / mean(dta$weighting_factor_orig)
# fivenum(dta$weighting_factor)


# observations with weight = 0 would have no influence, but we can't take their
# logarithm -> remove them from data set
dta <- dta[weighting_factor > 0]

# prepare data appropriately to be used in Poisson model
dta_grouped <- dta[, group_id := .GRP,
                   keyby = .(West_East, #region,
                             c_age,
                             West_East_c_age,
                             syear)][, .(share,
                                         West_East,
                                         # region,
                                         c_age,
                                         West_East_c_age,
                                         syear,
                                         weighting_factor,
                                         group_id)]

# Construct data set containing (unweighted) counts and weighted counts
step_size <- 0.01
n_bins <- 1/ step_size

dta_est <- construct_dataset(dta_grouped, response = list("share"),
                             step_size = step_size)
dta_est_weight <- construct_dataset(dta_grouped, response = list("share"),
                             step_size = step_size, weighted = TRUE)
dta_est[, weighted_counts := dta_est_weight$counts]

# dta_est[weighted_counts == 0 & counts != 0]

# We have to specify the argument weights in gam and use an offset in the formula
# to include the SOEP-weights properly; For each bin, the weight to pass to gam
# is the weighted count of observations in this bin, divided by the usual count
# (i.e., the number) of observations in this bin; the offset that has to be used
# is -log of the respective weight
dta_est[, gam_weights := ifelse(counts != 0, weighted_counts / counts, 1)] # If counts == 0, the value of weight can actually be arbitrary, but its log has to exist
dta_est[, gam_offsets := - log(gam_weights)]

# Check whether weights are approximately normalized: mean(dta_est$gam_weights)

# # Sanity check: number of covariate combinations
# max(dta_est$group_id)

# weights for discrete points are 1, equal equidistant stepsize for all covariate
# combinations (histograms) in continuous component
dta_est[, Delta := ifelse(dta_est$share %in% 0:1, 1, step_size)]

# In East Germany, observations start only from 1991 (1984 in West). However, when
# constructing a (penalized) B-spline basis for the interaction of syear with East,
# the whole range of syear is used, starting from 1984. Thus, for East Germany,
# we need a modified syear-variable, which contains the values of syear for East
# Germany and some value larger than 1991 for West Germany. (Note that since
# West Germany is the reference, there is no effect estimated at all; Furthermore,
# for interaction with categorical covariates, each category has a separate (usually
# identical, which is however not reasonable in our case) design matrix)
dta_est$syear_ <- ifelse(dta_est$West_East == "East", dta_est$syear, mean(dta_est$syear))

# For the interaction of syear with West_East and c_age, we need two separate
# categorical variables due to the different range of years in East. Note that
# for West, we no effect should be estimated for the reference West_3 (no minor
# children), but also no effect should be estimated for East Germany, which gets
# a separate variable. Thus, we combine both to West_3_East and use it as
# reference, i.e., no effect is estimated for both.
dta_est$West_c_age <- factor(ifelse(dta_est$West_East_c_age %in% c("West_3", paste0("East_", 1:3)),
                                    "West_3_East", as.character(dta_est$West_East_c_age)),
                             levels = c("West_3_East", paste0("West_", 2:1)), ordered = TRUE)
dta_est$East_c_age <- factor(ifelse(grepl("East", dta_est$West_East_c_age),
                                    as.character(dta_est$West_East_c_age), "West"),
                             levels = c("West", paste0("East_", 3:1)), ordered = TRUE)

saveRDS(dta_est, "./experiments/Eva_SOEP/Objects/dta_est.rds")

# ################################################################################
# ################### Inspect counts per covariate combination ###################
# ################################################################################
# # Could intercepts be modeled smoothly in the covariate basis?
# counts_per_hist <- sapply(seq_len(max(dta_est$group_id)),
#                           function(i) sum(dta_grouped[group_id == i, weighting_factor]))
#
# ylim <- c(0, max(counts_per_hist))
# legend_text <- rep(c("West", "East"), each = 3)
# legend_text_2 <- paste0("c_age ", rep(c("0-6", "7-18", "other"), 2))
# xaxt <- c(rep("n", 2),  "s", rep("n", 3))
# yaxt <- rep(c("s", "n"), each = 3)
# up_down <- 4.1
# mar <- list(c(0, 4.1, up_down, 0), c(0, 4.1, 0, 0), c(up_down, 4.1, 0, 0),
#             c(0, 0, up_down, 4.1), c(0, 0, 0, 4.1), c(up_down, 0, 0, 4.1))
# path <- "./experiments/Eva_SOEP/Images/"
#
# pdf(paste0(path, "counts_per_covariate_combination.pdf"), width = 10, height = 7)
# layout(matrix(1:6, nrow = 3), heights = c(1, 0.78, 1))
# for (i in 1:6) {
#   par(mar = mar[[i]])
#   years <- ifelse(i <= 3, 1984, 1991):2016
#   index <- if (i <= 3) {
#     (i - 1) * 33 + 1:33
#     } else {
#       3 * 33 + (i - 4) * 26 + 1:26
#     }
#   plot(years, counts_per_hist[index], ylab = "", xlab = "", ylim = ylim,
#        xaxt = xaxt[i], yaxt = yaxt[i])
#        # main = paste0("Observaion numbers for ", names(counts_per_hist)[i]))
#   abline(v = c(2000, 2010), col = "grey", lty = 3)
#   legend(x = "topleft", legend = c(legend_text[i], legend_text_2[i]),
#          bty = "n", cex = 1)
# }
# axis(1, at = seq(1995, 2015, by = 5))
# mtext(text = "observation numbers per covariate combination", side = 2, line = 36.5,
#       at = ylim[2] + ylim[2] /2) # 3.2
# mtext(text = "year", side = 1, line = 3, at = 1990)
# dev.off()
#
# # With weights included, there are only few less severe jumps in count levels
# # and it looks like the intercepts could indeed be modeled smoothly in the
# # covariate basis. However, the response data given to the model is the unweighted
# # one, which is not smooth (see below), and the weights are only passed via the
# # weights argument and the offset. I'm not sure, whether smoothness of the weighted
# # counts is sufficient to model the intercepts smoothly in the current model
# # specification. Since computation time is not that bad, I would continue
# # estimating one intercept per covariate combination.
# ### Old (without including weights:)
# # Not smooth, there are severe jumps in count levels in 2000 and 2010 (the latter
# # most pronounced for couples without minor children); Comparing this with the
# # overall development of the SOEP sample sizes (http://companion.soep.de/Target%20Population%20and%20Samples/Development%20of%20Sample%20Sizes.html)
# # one sees that in 2000, there was a refreshment of the sample, while in 2010,
# # there was a new sample group added, namely, Low-Income Families.

################################################################################
################################ Estimate model ################################
################################################################################

# For East Germany, we need to adjust the covariate basis for the interaction with
# syear, such that each B-spline function (before applying the sum-to-zero
# constraint) has at least one observation in its support. However, we want the
# remaining basis functions to coincide with the basis functions for West Germany.
# For this purpose, consider the knots generated for the whole range of syear
# (starting from 1984):
knots_east <- smooth.construct(ti(syear, share, bs = c("ps", "md"),
                                  m = list(c(2,2), c(2,2)),
                                  k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE),
                               data = dta_est, knots = NULL)$margin[[1]]$knots
# The first 5 knots are smaller than 1991:
which(knots_east < 1991)
# Since the support of each basis function (for second order B-splines) spans 5
# knots, this means that there are no observations in the support of the first
# basis function. Removing the smallest knot corresponds to removing the first
# basis function.
knots_east <- knots_east[-1]

# !!! Warning !!!
# The following call of gam() is time-consuming! Running it on one of HU's
# servers (768 GB memory; 12 physical cores at 3.4 GHz) takes 3-3.5 hours.
# The resulting model object was thus saved and can be loaded via:
model_soep <- readRDS("./experiments/Eva_SOEP/Objects/model_all_covariates_but_region.rds")

# Notes on model specification:
# - bs = "md" uses our self-written smoother for mixed reference measure. It
# already includes the integrate-to-zero constraint and thus no additional
# centering shall be applied -> use ti(..., mc = FALSE), since s() and te() don't
# have an option to turn centering constraints off!
# - np = TRUE (default) implements a reparametrization, which is not desired in
# our case -> use ti(..., np = FALSE)
# - offset(log(Delta)) adds the necessary additive term in the predictor that
# includes binwidths/dirac weights into the estimation
model_soep <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                              mc = FALSE, np = FALSE) # default: no penalty for discrete component
                  + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                       mc = FALSE, np = FALSE, by = West_East)
                  + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                       mc = FALSE, np = FALSE, by = c_age)
                  + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                       mc = FALSE, np = FALSE, by = West_East_c_age)
                  + ti(syear, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE)
                  + ti(syear_, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE, by = West_East)
                  + ti(syear, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE, by = c_age)
                  + ti(syear, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE, by = West_c_age)
                  + ti(syear_, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE, by = East_c_age)
                  + as.factor(group_id) - 1 # no scalar intercept, but one intercept per covariate combination
                  + offset(log(Delta) + gam_offsets),
                  data = dta_est, weights = dta_est$gam_weights,
                  method = "REML", family = poisson())

# Missing compared to first paper:
# region effect (centered around respective West_East-effect)
# However, result in first paper was that regions vary only very little... Including
# them would lead to 552 (instead of now 177) intercepts per covariate combination
# that need to be estimated, again increasing computation time...

saveRDS(model_soep, "./experiments/Eva_SOEP/Objects/model_all_covariates_but_region.rds")

################################################################################
############### Get predictions on clr- and retransformed level ################
################################################################################

X <- model.matrix(model_soep)
# remove intercepts per covariate combination (not of interest for estimated densities)
intercepts <- which(grepl("group_id", colnames(X)))
X <- X[, -intercepts]
theta_hat <- model_soep$coefficients[-intercepts]
f_hat_clr <- X %*% theta_hat

share <- unique(dta_est$share)
f_hat_clr <- matrix(f_hat_clr, nrow = length(share), ncol = length(intercepts))

# For our plots we use matrices with 33 columns for each value of West_East and
# c_age, i.e., we insert columns corresponding to the years 1984-1990 for East
# Germany with NAs
source("./experiments/Eva_SOEP/help_functions.R")
pred_clr_j <- fill_missing_years(f_hat_clr) # Add NA columns for the
# Due to reference coding, order of c_age is 3, 2, 1 (in blocks of 33 columns,
# where each column corresponds to one year); For plots, we want 1, 2, 3
pred_clr_j <- pred_clr_j[, c(67:99, 34:66, 1:33, 166:198, 133:165, 100:132)]
saveRDS(pred_clr_j, "./experiments/Eva_SOEP/Objects/pred_clr_j.rds")

# Transform to Bayes-level
pred_j <- clr_inv(pred_clr_j)
# Problem (after including weights): t(pred_j[, 1:2]) contains a lot of NaNs;
# Reason: On clr-level t(pred_clr_j[, 1:2]) contains a lot of large values, whose
# exponential is Inf in R, see t(exp(pred_clr_j[, 1:2])), furthermore, the
# denominator int exp(pred_clr_j) dmu is Inf and Inf/Inf results in NaN...
saveRDS(pred_j, "./experiments/Eva_SOEP/Objects/pred_j.rds")

################################################################################
###################### Save "histograms" for later plots #######################
################################################################################
# For each covariate combination: Compute underlying histogram on relative (i.e.,
# density) scale
# without weights: nrow(dta_grouped[group_id == g]) instead of sum(dta_grouped[group_id == g, weighting_factor])
zeros_dens <- sapply(unique(dta_est$group_id),
                     function(g) dta_est[group_id == g & share == 0]$weighted_counts / sum(dta_grouped[group_id == g, weighting_factor]))
dual_earner_dens <- sapply(unique(dta_est$group_id),
                           function(g) {
                             counts_dens <- dta_est[group_id == g & !(share %in% 0:1)]$weighted_counts
                             sum(counts_dens) / sum(dta_grouped[group_id == g, weighting_factor]) * counts_dens / (step_size * sum(counts_dens))
                           })
ones_dens <- sapply(unique(dta_est$group_id),
                    function(g) dta_est[group_id == g & share == 1]$weighted_counts / sum(dta_grouped[group_id == g, weighting_factor]))

# Due to reference coding, order of c_age is 3, 2, 1 (in blocks of 33 columns,
# where each column corresponds to one year); For plots, we want 1, 2, 3
rearrange <- function(x) {
  order_with_missings <- c(67:99, 34:66, 1:33, 166:198, 133:165, 100:132)
  order_without_missings <- c(67:99, 34:66, 1:33, 152:177, 126:151, 100:125)
  if (length(x) == 6 * 33) {
    return(x[order_with_missings])
  } else if (length(x) == 3 * 33 + 3 * 26) {
    return(x[order_without_missings])
  } else {
    stop("x has wrong length")
  }
}
zeros_dens <- rearrange(zeros_dens)
dual_earner_dens <- t(apply(dual_earner_dens, 1, rearrange))
ones_dens <- rearrange(ones_dens)

hist_data <- list(zeros_dens = zeros_dens, dual_earner_dens = dual_earner_dens,
                  ones_dens = ones_dens, share = share)

# # Check:
# sapply(seq_along(hist_data$zeros_dens),
#        function(i) hist_data$zeros_dens[i] + 0.01 * sum(hist_data$dual_earner_dens[, i]) + hist_data$ones_dens[i])

saveRDS(hist_data, "./experiments/Eva_SOEP/Objects/hist_data.rds")
