current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/../.."))

# Package to estimate GAMs (like our Poisson model)
library(mgcv)

model_soep <- readRDS("./experiments/Eva_SOEP/Objects/model_all_covariates_but_region.rds")
path <- "./experiments/Eva_SOEP/Images/"

pdf(paste0(path, "compare_smoothing_parameters.pdf"), width = 12, height = 8)
par(mar = c(17, 4, 4, 2) + 0.1, mfrow = c(1, 2))
plot(1:length(model_soep$sp), model_soep$sp, xaxt = "n", xlab = "",
     ylab = "smoothing parameter", pch = 4)
axis(1, at = 1:length(model_soep$sp), labels = names(model_soep$sp), las = 2)

plot(1:length(model_soep$sp), model_soep$sp, xaxt = "n", xlab = "",
     ylim = c(0, 1), ylab = "smoothing parameter", pch = 4)
axis(1, at = 1:length(model_soep$sp), labels = names(model_soep$sp), las = 2)
dev.off()

##### Comparison of old models (saved in ./experiments/Eva_SOEP/Reparametrization_500/Images/)
model_soep_500 <- readRDS("./experiments/Eva_SOEP/Reparametrization_500/Objects/model_all_covariates_but_region.rds")
model_soep <- readRDS("./experiments/Eva_SOEP/No_reparametrization/Objects/model_all_covariates_but_region.rds")
model_soep_wo_weight <- readRDS("./experiments/Eva_SOEP/Without_weights/Objects/model_all_covariates_but_region.rds")

path <- "./experiments/Eva_SOEP/Reparametrization_500/Images/"

pdf(paste0(path, "compare_smoothing_parameters.pdf"), width = 12, height = 8)
par(mar = c(17, 4, 4, 2) + 0.1, mfrow = c(1, 2))
plot(1:length(model_soep$sp), model_soep_500$sp, xaxt = "n", xlab = "",
     ylim = range(model_soep_500$sp, model_soep_wo_weight$sp, model_soep$sp),
     ylab = "smoothing parameter", pch = 4)
points(1:length(model_soep$sp), model_soep$sp, col = 2, pch = 4)
points(1:length(model_soep$sp), model_soep_wo_weight$sp, col = 3, pch = 4)
axis(1, at = 1:length(model_soep$sp), labels = names(model_soep$sp), las = 2)
# abline(h = sort(c(model_soep_500$sp, model_soep_wo_weight$sp, model_soep$sp), decreasing = TRUE)[10],
#        col = "grey", lty = 3)
legend(1, 270000, legend = c("Rescaled weights", "Original weights", "No weights"),
       pch = 4, col = 1:3, bty = "n", xpd = TRUE, horiz = TRUE)

# # Log-scale
# plot(1:length(model_soep$sp), log10(model_soep_500$sp), xaxt = "n", xlab = "",
#      ylim = log10(range(model_soep_500$sp, model_soep_wo_weight$sp, model_soep$sp)),
#      ylab = expression(paste(log[10], "(smoothing parameter)")), pch = 4)
# points(1:length(model_soep$sp), log10(model_soep$sp), col = 2, pch = 4)
# points(1:length(model_soep$sp), log10(model_soep_wo_weight$sp), col = 3, pch = 4)
# axis(1, at = 1:length(model_soep$sp), labels = names(model_soep$sp), las = 2)

plot(1:length(model_soep$sp), model_soep_500$sp, xaxt = "n", xlab = "",
     ylim = c(0, 1), # range(sort(c(model_soep_500$sp, model_soep_wo_weight$sp, model_soep$sp), decreasing = TRUE)[-(1:9)]),
     ylab = "smoothing parameter", pch = 4)
points(1:length(model_soep$sp), model_soep$sp, col = 2, pch = 4)
points(1:length(model_soep$sp), model_soep_wo_weight$sp, col = 3, pch = 4)
axis(1, at = 1:length(model_soep$sp), labels = names(model_soep$sp), las = 2)
# abline(h = sort(c(model_soep_500$sp, model_soep_wo_weight$sp, model_soep$sp), decreasing = TRUE)[10],
#        col = "grey", lty = 3)
dev.off()
