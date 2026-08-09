setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# devtools::install_github("Eva2703/DensityRegression")
library(DensityRegression)
library(data.table)

scipen <- getOption("scipen") # save current scipen setting for resetting it later
options(scipen = 999) # preventing scientific notation (in particular when saving objects) for large sample size

# source("help_functions.R")
path <- "./Objects/"

source("evaluate_sim_help_functions.R")

# model_orig <- readRDS("../../Eva_SOEP/Objects/model_soep_dens_Reg_mc_syear__TRUE.rds")
# model_orig <- readRDS("../../Eva_SOEP/Objects/model_soep_dens_Reg.rds")

# Check, whether in any iteration, simulation model object was saved (which is only
# done, if the design matrices in the simulation model do not match
# the oned in the true model)
length(list.files(path = path, pattern = "^model_sim_.*.rds$", full.names = TRUE))
# 0 -> all models use correct design matrix

N <- c(50000, 150000, 500000)
G <- c(50, 100, 200) # c(25, 50, 100, 200)

################################################################################
############################### Evaluate relMSE ################################
################################################################################
### predictions
relMSE_predictions <- list()
for (i in seq_along(G)) {
  relMSE_predictions[[i]] <- list()
  for (j in seq_along(N)) {
    relMSE_predictions[[i]][[j]] <- list.files(path = path,
                                               pattern = paste0("relMSE_N_predictions_.*N_",
                                                                N[j], "_G_", G[i], ".rds"),
                                               full.names = TRUE)
    relMSE_predictions[[i]][[j]] <- sapply(relMSE_predictions[[i]][[j]], readRDS)
    names(relMSE_predictions[[i]][[j]]) <- stringr::str_remove_all(names(relMSE_predictions[[i]][[j]]), 
                                                                   "(./Objects/relMSE_N_predictions_|.rds)")
  }
}

# # Check completeness
# sapply(relMSE_predictions, lengths)
# done <- list()
# for (i in seq_along(G)) {
#   done[[i]] <- list()
#   for (j in seq_along(N)) {
#     done[[i]][[j]] <- list.files(path = path, pattern = paste0("relMSE_N_predictions_.*N_",
#                                                                 N[j], "_G_", G[i], ".rds"),
#                                  full.names = TRUE)
#     done[[i]][[j]] <- sort(as.numeric(stringr::str_extract(done[[i]][[j]], "[0-9]+")))
#   }
# }
# lapply(done, function(d) lapply(d, function(d_) which(!(1:200 %in% d_))))

# # par(mfrow = c(length(G), length(N)))
# ylims <- lapply(seq_along(N), 
#                 function(j) c(0, max(unlist(sapply(seq_along(G), 
#                                                    function(i) relMSE_predictions[[i]][[j]])))))

# ylims[[1]] <- ylims[[2]]

# # Outlier detection:
# # Which iteratations yield the largest relMSE, respectively? -> Noch mal anschauen, wenn alles durchgelaufen
# lapply(seq_along(relMSE_predictions), 
#        function(i) sapply(seq_along(relMSE_predictions[[i]]), 
#                           function(j) which(relMSE_predictions[[i]][[j]] == 
#                                               max(relMSE_predictions[[i]][[j]]))))
# # How large do the relMSEs get?
# lapply(seq_along(relMSE_predictions), 
#        function(i) lengths(sapply(seq_along(relMSE_predictions[[i]]), 
#                           function(j) which(relMSE_predictions[[i]][[j]] > 1.5))))
# # # Largest for smallest N (50000)
# # lapply(seq_along(relMSE_predictions), 
# #        function(i) sapply(seq_along(relMSE_predictions[[i]]), 
# #                           function(j) sort(unlist(relMSE_predictions[[i]][[j]]), decreasing = TRUE)[1:15]))
# sapply(seq_along(relMSE_predictions),
#        function(i) sort(relMSE_predictions[[i]][[1]], decreasing = TRUE)[1:20])

# # xlims <- lapply(relMSE_predictions, function(relMSE) range(unlist(relMSE)))
# xlims <- lapply(relMSE_predictions, function(relMSE) c(0, min(max(unlist(relMSE)), 1.5)))
# params <- get_params_for_matrix_plot(n_cols = length(G), n_rows = length(N), byrow = FALSE)
# def.par <- par(no.readonly = TRUE)
# 
# pdf("./Images/relMSE_pred_matrix.pdf")
# layout(matrix(1:(length(N) * length(G)), nrow = length(N)), 
#        heights = c(1, rep(0.74, length(N) - 2), 1), 
#        widths = c(1, rep(0.85, length(G) - 2), 1))
# sapply(seq_along(relMSE_predictions), 
#        function(i) sapply(seq_along(relMSE_predictions[[i]]), 
#                           function(j) { 
#                             param_ind <- (i - 1) * length(N) + j
#                             if (j == 1) {
#                               main <- paste0("G = ", G[i])
#                             } else {
#                               main <- ""
#                             }
#                             par(mar = params$mar[[param_ind]])
#                             boxplot(relMSE_predictions[[i]][[j]], 
#                                     main = main, horizontal = TRUE,
#                                     # main = paste0("G = ", G[i], ", N = ", N[j]),
#                                     ylim = xlims[[i]], # ylim = ylims[[j]],
#                                     xaxt = params$xaxt[param_ind], yaxt = params$yaxt[param_ind]
#                                     )
#                             if (i == length(G)) {
#                               mtext(paste0("N = ", N[j]), side = 4, line = 1, 
#                                     las = 1, cex = 0.8)
#                             }
#                             }))
# dev.off()
# 
# lapply(seq_along(relMSE_predictions),
#        function(i) sapply(seq_along(relMSE_predictions[[i]]),
#                           function(j) fivenum(unlist(relMSE_predictions[[i]][[j]]))))
# params <- get_params_for_matrix_plot(n_cols = length(G), n_rows = 1, byrow = FALSE,
#                                      up = 3, le_ri = 5.5)
# # def.par <- par(no.readonly = TRUE)
# 
# # pdf("./Images/relMSE_first_version.pdf")
# layout(matrix(1:(length(G) * 1), ncol = length(G)), 
#        # heights = c(1, rep(0.74, length(N) - 2), 1), 
#        widths = c(1, rep(0.85, length(G) - 2), 1))
# sapply(seq_along(relMSE_predictions), 
#        function(i) {
#          main <- paste0("G = ", G[i])
#          par(mar = params$mar[[i]])
#          boxplot(relMSE_predictions[[i]][length(relMSE_predictions[[i]]):1], 
#                  # main = main, 
#                  horizontal = TRUE,
#                  # main = paste0("G = ", G[i], ", N = ", N[j]),
#                  ylim = xlims[[i]], # ylim = ylims[[j]],
#                  xaxt = params$xaxt[i], yaxt = "n"
#          )
#          mtext(text = main, side = 3, line = 1)
#          if (i == length(G)) {
#            axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
#          }
#        }
#        )

### partial effects
relMSE_effects <- list()
for (i in seq_along(G)) {
  relMSE_effects[[i]] <- list()
  for (j in seq_along(N)) {
    relMSE_effects[[i]][[j]] <- list.files(path = path, 
                                           pattern = paste0("relMSE_effects_N_.*N_", 
                                                            N[j], "_G_", G[i], ".rds"),
                                           full.names = TRUE)
    relMSE_effects[[i]][[j]] <- sapply(relMSE_effects[[i]][[j]], readRDS)
    colnames(relMSE_effects[[i]][[j]]) <- stringr::str_remove_all(colnames(relMSE_effects[[i]][[j]]), 
                                                                      "(./Objects/relMSE_effects_N_|.rds)")
  }
}

# # simultaneous effects correspond to first 8 rows, respectively
# # rownames(relMSE_effects[[i]][[1]])
# N_effects <- 8 # nrow(relMSE_effects[[1]][[1]])
# 
# library(latex2exp)
# effect_names <- c(TeX("$\\hat{\\beta}_0$"), TeX("$\\hat{\\beta}_{West\\_East}$"), 
#                   TeX("$\\hat{\\beta}_{c\\_age}$"), TeX("$\\hat{\\beta}_{c\\_age, West\\_East}$"),
#                   TeX("$\\hat{g}(year)$"), TeX("$\\hat{g}_{West\\_East}(year)$"),
#                   TeX("$\\hat{g}_{c\\_age}(year)$"), TeX("$\\hat{g}_{c\\_age, West\\_East}(year)$"),
#                   TeX("$\\hat{\\beta}_{2}$"), TeX("$\\hat{\\beta}_{1}$"),
#                   TeX("$\\hat{\\beta}_{2, East}$"), TeX("$\\hat{\\beta}_{1, East}$"),
#                   TeX("$\\hat{g}_{2}(year)$"), TeX("$\\hat{g}_{1}(year)$"),
#                   TeX("$\\hat{g}_{2, East}(year)$"), TeX("$\\hat{g}_{1, East}(year)$"))
# params <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_effects,
#                                      byrow = FALSE, up = 3, le_ri = 6.5)
# def.par <- par(no.readonly = TRUE)

# test <- relMSE_effects[[1]]
# lapply(seq_along(test), function(i) test[[i]][1,])
# test_ <- lapply(seq_along(G), 
#                 function(i) lapply(seq_len(N_effects), 
#                                    function(k) lapply(seq_along(N), 
#                                                       function(j) relMSE_effects[[i]][[j]][k,])))

# pdf("./Images/relMSE_effects.pdf", height = 20, width = 7)
# layout(matrix(1:(length(G) * N_effects), ncol = length(G)),
#        heights = c(1, rep(0.54, N_effects - 2), 1),
#        widths = c(1, rep(0.85, length(G) - 2), 1))
# sapply(seq_along(relMSE_effects),
#        function(i)
#          sapply(seq_len(N_effects),
#                 function(k) {
#                   param_ind <- (i - 1) * N_effects + k
#                   if (k == 1) {
#                     main <- paste0("G = ", G[i])
#                   } else {
#                     main <- ""
#                   }
#                   par(mar = params$mar[[param_ind]])
#                   boxplot(lapply(seq_along(N),
#                                  function(j) relMSE_effects[[i]][[j]][k,])[length(N):1],
#                           # main = main,
#                           horizontal = TRUE,
#                           # main = paste0("G = ", G[i], ", N = ", N[j]),
#                           # ylim = xlims[[i]], # ylim = ylims[[j]],
#                           ylim = c(0, 3),
#                           xaxt = params$xaxt[param_ind], yaxt = "n"
#                   )
#                   mtext(text = main, side = 3, line = 1)
#                   if (i == length(G)) {
#                     axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
#                   }
#                   if (i == 1) {
#                     mtext(text = effect_names[k], # rownames(relMSE_effects[[i]][[1]])[k],
#                           side = 2, las = 1, line = 1)
#                   }
#                 }
#          ))
# # dev.off()

### predictions and partial effects combined

relMSE <- lapply(seq_along(G), 
                 function(i) lapply(seq_along(N), 
                                    function(j) rbind("pred" = relMSE_predictions[[i]][[j]], relMSE_effects[[i]][[j]])))

# lapply(seq_along(relMSE),
#        function(i) lapply(seq_len(N_effects),
#                           function(k) sapply(seq_along(N), 
#                                              function(j) fivenum((relMSE[[i]][[j]][k,])))))

# 1st row = prediction, rows 2-9: simultaneous effects (rows 10-16: pointwise effect, where different)
# rownames(relMSE[[i]][[1]])
# N_effects <- 9 # nrow(relMSE[[1]][[1]])

### Only main effects (main paper)
main_effects <- c(1:6)[-5]
N_main_effects <- length(main_effects)

library(latex2exp)
effect_names <- c(TeX("$\\hat{\\f}$"), TeX("$\\hat{\\beta}_0$"), 
                  TeX("$\\hat{\\beta}_{West\\_East}$"), TeX("$\\hat{\\beta}_{c\\_age}$"), 
                  TeX("$\\hat{\\beta}_{c\\_age, West\\_East}$"),
                  TeX("$\\hat{g}(year)$"), TeX("$\\hat{g}_{West\\_East}(year)$"),
                  TeX("$\\hat{g}_{c\\_age}(year)$"), TeX("$\\hat{g}_{c\\_age, West\\_East}(year)$"),
                  TeX("$\\hat{\\beta}_{2}$"), TeX("$\\hat{\\beta}_{1}$"),
                  TeX("$\\hat{\\beta}_{2, East}$"), TeX("$\\hat{\\beta}_{1, East}$"),
                  TeX("$\\hat{g}_{2}(year)$"), TeX("$\\hat{g}_{1}(year)$"),
                  TeX("$\\hat{g}_{2, East}(year)$"), TeX("$\\hat{g}_{1, East}(year)$"))


lapply(seq_along(relMSE),
       function(i) lapply(seq_len(N_main_effects),
                          function(k) sapply(seq_along(N), 
                                             function(j) fivenum((relMSE[[i]][[j]][main_effects[k],])))))

lapply(seq_along(relMSE),
       function(i) lapply(seq_len(N_main_effects),
                          function(k) sort((relMSE[[i]][[1]][main_effects[k],]), decreasing = TRUE)[1:20]))

ymax_main <- 3.8
# To add number of outliers, we cannot plot the whole range
# sort(relMSE[[3]][[1]][4,], decreasing = TRUE)[1:20]
ycut_main <- 3.35

params_main <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_main_effects,
                                          byrow = FALSE, up = 3, le_ri = 8.5)
pdf("./Images/relMSE_main.pdf", height = 4, width = 6.5)
layout(matrix(1:(length(G) * N_main_effects), ncol = length(G)),
       heights = c(1, rep(0.47, N_main_effects - 2), 1),
       widths = c(1, rep(0.56, length(G) - 2), 1))
sapply(seq_along(relMSE),
       function(i)
         sapply(seq_len(N_main_effects),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  relMSEs <- lapply(seq_along(N),
                                    function(j) relMSE[[i]][[j]][main_effects[k],])
                  out_of_lim <- lapply(relMSEs, function(r) length(which(r > ycut_main)))
                  relMSEs_plot <- lapply(relMSEs,
                                         function(r) r[which(r <= ycut_main)])
                  par(mar = params_main$mar[[param_ind]])
                  boxplot(relMSEs_plot[length(N):1],
                          # main = main,
                          horizontal = TRUE, lwd = 0.6,
                          # main = paste0("G = ", G[i], ", N = ", N[j]),
                          # ylim = xlims[[i]], # ylim = ylims[[j]],
                          ylim = c(0, ymax_main),
                          xaxt = params_main$xaxt[param_ind], yaxt = "n"
                  )
                  abline(v = ycut_main, lty = 1, lwd = 0.6)
                  sapply(seq_along(out_of_lim), function(n) {
                    if (out_of_lim[[n]] > 0) {
                      text(x = ymax_main * 0.95, y = (length(out_of_lim):1)[n],
                           labels = paste0("+", out_of_lim[n])) # , col = "red")
                    }
                  })
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[main_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  # if (k == N_main_effects) {
                  #   segments(x0 = ycut_main, y0 = 0.5, x1 = ycut_main, y1 = -1.1,
                  #            lwd = 0.6, xpd = TRUE)
                  # }
                }
         ))
at_y <- 3 + 3 * (N_main_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * ymax_main * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$relMSE(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

### All effects (appendix)
all_effects <- 1:9
N_all_effects <- length(all_effects)


# lapply(seq_along(relMSE),
#        function(i) lapply(seq_len(N_all_effects),
#                           function(k) sapply(seq_along(N), 
#                                              function(j) fivenum((relMSE[[i]][[j]][all_effects[k],])))))
# 
# lapply(seq_along(relMSE),
#        function(i) lapply(seq_len(N_all_effects),
#                           function(k) sort((relMSE[[i]][[1]][all_effects[k],]), decreasing = TRUE)[1:20]))

ymax_all <- 15
# sort(relMSE[[3]][[1]][8,], decreasing = TRUE)[1:20]
ycut_all <- 13

params_all <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_all_effects,
                                     byrow = FALSE, up = 3, le_ri = 13.5)
pdf("./Images/relMSE_all.pdf", height = 6.3, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(relMSE),
       function(i)
         sapply(seq_len(N_all_effects),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  relMSEs <- lapply(seq_along(N),
                                    function(j) relMSE[[i]][[j]][all_effects[k],])
                  out_of_lim <- lapply(relMSEs, function(r) length(which(r > ycut_all)))
                  relMSEs_plot <- lapply(relMSEs,
                                         function(r) r[which(r <= ycut_all)])
                  par(mar = params_all$mar[[param_ind]])
                  boxplot(relMSEs_plot[length(N):1],
                          # main = main,
                          horizontal = TRUE, lwd = 0.6,
                          # main = paste0("G = ", G[i], ", N = ", N[j]),
                          # ylim = xlims[[i]], # ylim = ylims[[j]],
                          ylim = c(0, ymax_all),
                          xaxt = "n", # params_all$xaxt[param_ind], 
                          yaxt = "n"
                  )
                  abline(v = ycut_all, lty = 1, lwd = 0.6)
                  sapply(seq_along(out_of_lim), function(n) {
                    if (out_of_lim[[n]] > 0) {
                      text(x = ymax_all * 0.95, y = (length(out_of_lim):1)[n],
                           labels = paste0("+", out_of_lim[n])) # , col = "red")
                    }
                  })
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[all_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all$xaxt[param_ind] == "s") {
                    axis(side = 1, at = c(0, 5, 10))
                  }
                  # if (k == N_all_effects) {
                  #   segments(x0 = ycut_main, y0 = 0.5, x1 = ycut_main, y1 = -1.1,
                  #            lwd = 0.6, xpd = TRUE)
                  # }
                }
         ))
at_y <- 3 + 3 * (N_all_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * ymax_all * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$relMSE(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

### MSE with denominator
MSE_predictions <- list()
for (i in seq_along(G)) {
  MSE_predictions[[i]] <- list()
  for (j in seq_along(N)) {
    MSE_predictions[[i]][[j]] <- list.files(path = path,
                                               pattern = paste0("^MSE_N_predictions_.*N_",
                                                                N[j], "_G_", G[i], ".rds"),
                                               full.names = TRUE)
    MSE_predictions[[i]][[j]] <- sapply(MSE_predictions[[i]][[j]], readRDS)
    names(MSE_predictions[[i]][[j]]) <- stringr::str_remove_all(names(MSE_predictions[[i]][[j]]), 
                                                                   "(./Objects/MSE_N_predictions_|.rds)")
  }
}

MSE_effects <- list()
for (i in seq_along(G)) {
  MSE_effects[[i]] <- list()
  for (j in seq_along(N)) {
    MSE_effects[[i]][[j]] <- list.files(path = path, 
                                        pattern = paste0("^MSE_effects_N_.*N_", 
                                                         N[j], "_G_", G[i], ".rds"),
                                        full.names = TRUE)
    MSE_effects[[i]][[j]] <- sapply(MSE_effects[[i]][[j]], readRDS)
    colnames(MSE_effects[[i]][[j]]) <- stringr::str_remove_all(colnames(MSE_effects[[i]][[j]]), 
                                                               "(./Objects/MSE_effects_N_|.rds)")
  }
}
MSE <- lapply(seq_along(G), 
                 function(i) lapply(seq_along(N), 
                                    function(j) rbind("pred" = MSE_predictions[[i]][[j]], MSE_effects[[i]][[j]])))

lapply(seq_along(MSE),
       function(i) lapply(seq_len(N_main_effects),
                          function(k) sapply(seq_along(N), 
                                             function(j) fivenum((MSE[[i]][[j]][k,])))))

denominator_effects_L <- readRDS(paste0(path, "denominator_effects_relMSE_L.rds"))
average_categories <- function(x) {
  c(x[1:2], mean(x[3:4]), mean(x[5:6]), x[7:8], mean(x[9:10]), mean(x[11:12]))
}
denominator_effects_L_agg <- average_categories(denominator_effects_L)
denominator_L <- readRDS(paste0(path, "denominator_relMSE_L.rds"))
denominator <- c(denominator_L, denominator_effects_L_agg)

# More insights on year-effect for c_age (which yields largest relMSEs):
norm_truth_effects <- readRDS(paste0(path, "norm_truth_effects.rds"))
fivenum(unlist(norm_truth_effects[9:10]))
mean(unlist(norm_truth_effects[9:10])) == denominator[8]

ymax_all_MSE <- 30
# sort(MSE[[2]][[1]][9,], decreasing = TRUE)[1:20]
ycut_all_MSE <- 26 # 25.8

params_all_MSE <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_all_effects,
                                     byrow = FALSE, up = 3, le_ri = 13.5)
pdf("./Images/MSE_all.pdf", height = 6.3, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(MSE),
       function(i)
         sapply(seq_len(N_all_effects),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  MSEs <- lapply(seq_along(N),
                                 function(j) MSE[[i]][[j]][all_effects[k],])
                  out_of_lim <- lapply(MSEs, function(r) length(which(r > ycut_all_MSE)))
                  MSEs_plot <- lapply(MSEs,
                                      function(r) r[which(r <= ycut_all_MSE)])
                  par(mar = params_all_MSE$mar[[param_ind]])
                  boxplot(MSEs_plot[length(N):1],
                          # main = main,
                          horizontal = TRUE, lwd = 0.6,
                          # main = paste0("G = ", G[i], ", N = ", N[j]),
                          # ylim = xlims[[i]], # ylim = ylims[[j]],
                          ylim = c(0, ymax_all_MSE),
                          xaxt = "n", # params_all_MSE$xaxt[param_ind], 
                          yaxt = "n"
                  )
                  abline(v = ycut_all_MSE, lty = 1, lwd = 0.6)
                  sapply(seq_along(out_of_lim), function(n) {
                    if (out_of_lim[[n]] > 0) {
                      text(x = ymax_all_MSE * 0.95, y = (length(out_of_lim):1)[n],
                           labels = paste0("+", out_of_lim[n])) # , col = "red")
                    }
                  })
                  # points(x = denominator[k], y = 1:3)
                  abline(v = denominator[k], col = "deepskyblue", lty = 2, lwd = 0.9)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[all_effects[k]], # rownames(MSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 25, by = 5))
                  }
                  # if (k == N_all_effects) {
                  #   segments(x0 = ycut_all_MSE, y0 = 0.5, x1 = ycut_all, y1 = -1.1,
                  #            lwd = 0.6, xpd = TRUE)
                  # }
                }
         ))
at_y <- 3 + 3 * (N_all_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * ymax_all_MSE * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$MSE(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

################################################################################
############################## Evaluate coverage ###############################
################################################################################
### Predictions
# Vc
coverage_predictions_Vc <- list()
empirical_coverage_predictions_Vc <- list()
for (i in seq_along(G)) {
  coverage_predictions_Vc[[i]] <- list()
  empirical_coverage_predictions_Vc[[i]] <- list()
  for (j in seq_along(N)) {
    coverage_predictions_Vc[[i]][[j]] <- list.files(path = path,
                                                    pattern = paste0("coverage_predictions_Vc_.*N_",
                                                                     N[j], "_G_", G[i], ".rds"),
                                                    full.names = TRUE)
    coverage_predictions_Vc[[i]][[j]] <- lapply(coverage_predictions_Vc[[i]][[j]], readRDS)
    empirical_coverage_predictions_Vc[[i]][[j]] <- apply(matrix(unlist(coverage_predictions_Vc[[i]][[j]]),
                                                                byrow = TRUE,
                                                                nrow = length(coverage_predictions_Vc[[i]][[j]])),
                                                         2, mean)
  }
}

# Plot coverage per conditional density (Vc)
params_pred <- get_params_for_matrix_plot(n_cols = length(G), n_rows = length(N),
                                          byrow = FALSE, up = 3, le_ri = 13.5)
hxpos <- c(rep(-0.7, 2), -0.9, rep(-0.7, 2))
pdf("./Images/coverage_Vc_pred.pdf", height = 8, width = 7.8) # width = 6.5
layout(matrix(1:(length(G) * length(N)), ncol = length(G)),
       heights = c(1, rep(0.78, length(N) - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_along(N),
                function(j) {
                  param_ind <- (i - 1) * length(N) + j
                  if (j == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_pred$mar[[param_ind]])
                  plot(empirical_coverage_predictions_Vc[[i]][[j]],
                       -(1:177), xlim = c(0, 1), xlab = "", ylab = "", ylim = c(-180.5, 3.5),
                       xaxt = "n", # params_main$xaxt[param_ind], 
                       yaxt =  "n", yaxs = "i", pch = 20)
                  lines(empirical_coverage_predictions_Vc[[i]][[j]][1:33], -(1:33), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vc[[i]][[j]][34:66], -(34:66), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vc[[i]][[j]][67:99], -(67:99), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vc[[i]][[j]][100:125], -(100:125), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vc[[i]][[j]][126:151], -(126:151), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vc[[i]][[j]][152:177], -(152:177), 
                        lwd = 0.6)
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  hlines <- -c(33.5, 66.5, 99.5, 99.5 + 26, 99.5 + 2 * 26)
                  abline(h = hlines, lwd = 0.6)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = -88.5, text = paste0("N = ", N[j]), 
                          las = 1, line = 1)
                  }
                  # if (i == 1) {
                  #   mtext(text = effect_names[main_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                  #         side = 2, las = 1, line = 1)
                  # }
                  if (params_pred$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                  if (params_pred$yaxt[param_ind] == "s") {
                    yaxis_at <- c(seq(2, 32, by = 5), 33 + seq(2, 32, by = 5), 
                                  66 + seq(2, 32, by = 5), 99 + seq(5, 25, by = 5),
                                  99 + 26 + seq(5, 25, by = 5), 99 + 2 * 26 + seq(5, 25, by = 5))
                    axis(side = 2, at = -yaxis_at,
                         labels = c(rep(1984:2016, 3), rep(1991:2016, 3))[yaxis_at],
                         las = 1)
                    mtext(text = rep(c("other", "7-18", "0-6"), 2), side = 2,
                          at = -c(17, 17 + 33, 17 + 66, 99 + 13.6, 
                                  99 + 13.6 + 26, 99 + 13.6 + 2 * 26),
                          line = 4, las = 1, cex = 0.75)
                    mtext(text = c("West", "East"), side = 2,
                          at = -c(17 + 33, 99 + 13.6 + 26),
                          line = 7.3, # las = 1, 
                          cex = 0.75)
                    sapply(seq_along(hlines), 
                           function(h) segments(x0 = 0, y0 = hlines[h], x1 = hxpos[h], 
                                                y1 = hlines[h], lwd = 0.6, xpd = TRUE))
                  }
                  if (i == 1 & j %in% 1:2) {
                    segments(x0 = 0, y0 = -180.5, x1 = -1.1, 
                             y1 = -180.5, xpd = TRUE)
                  }
                }
         ))
at_y <- 177 / 2 * 1.1 # 3 + 3 * (N_main_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{f}_{West\\_East, c\\_age, year}$"), at = at_y, line = 32.5, side = 2)
mtext(text = TeX("$empCR(\\hat{f}_{West\\_East, c\\_age, year})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# Vp
coverage_predictions_Vp <- list()
empirical_coverage_predictions_Vp <- list()
for (i in seq_along(G)) {
  coverage_predictions_Vp[[i]] <- list()
  empirical_coverage_predictions_Vp[[i]] <- list()
  for (j in seq_along(N)) {
    coverage_predictions_Vp[[i]][[j]] <- list.files(path = path,
                                                    pattern = paste0("coverage_predictions_Vp_.*N_",
                                                                     N[j], "_G_", G[i], ".rds"),
                                                    full.names = TRUE)
    coverage_predictions_Vp[[i]][[j]] <- lapply(coverage_predictions_Vp[[i]][[j]], readRDS)
    empirical_coverage_predictions_Vp[[i]][[j]] <- apply(matrix(unlist(coverage_predictions_Vp[[i]][[j]]),
                                                                byrow = TRUE,
                                                                nrow = length(coverage_predictions_Vp[[i]][[j]])),
                                                         2, mean)
  }
}

# Plot coverage per conditional density (Vp)
# params_pred <- get_params_for_matrix_plot(n_cols = length(G), n_rows = length(N),
#                                           byrow = FALSE, up = 3, le_ri = 13.5)
# hxpos <- c(rep(-0.7, 2), -0.9, rep(-0.7, 2))
pdf("./Images/coverage_Vp_pred.pdf", height = 8, width = 7.8) # width = 6.5
layout(matrix(1:(length(G) * length(N)), ncol = length(G)),
       heights = c(1, rep(0.78, length(N) - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_along(N),
                function(j) {
                  param_ind <- (i - 1) * length(N) + j
                  if (j == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_pred$mar[[param_ind]])
                  plot(empirical_coverage_predictions_Vp[[i]][[j]],
                       -(1:177), xlim = c(0, 1), xlab = "", ylab = "", ylim = c(-180.5, 3.5),
                       xaxt = "n", # params_main$xaxt[param_ind], 
                       yaxt =  "n", yaxs = "i", pch = 20)
                  lines(empirical_coverage_predictions_Vp[[i]][[j]][1:33], -(1:33), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vp[[i]][[j]][34:66], -(34:66), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vp[[i]][[j]][67:99], -(67:99), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vp[[i]][[j]][100:125], -(100:125), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vp[[i]][[j]][126:151], -(126:151), 
                        lwd = 0.6)
                  lines(empirical_coverage_predictions_Vp[[i]][[j]][152:177], -(152:177), 
                        lwd = 0.6)
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  hlines <- -c(33.5, 66.5, 99.5, 99.5 + 26, 99.5 + 2 * 26)
                  abline(h = hlines, lwd = 0.6)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = -88.5, text = paste0("N = ", N[j]), 
                          las = 1, line = 1)
                  }
                  # if (i == 1) {
                  #   mtext(text = effect_names[main_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                  #         side = 2, las = 1, line = 1)
                  # }
                  if (params_pred$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                  if (params_pred$yaxt[param_ind] == "s") {
                    yaxis_at <- c(seq(2, 32, by = 5), 33 + seq(2, 32, by = 5), 
                                 66 + seq(2, 32, by = 5), 99 + seq(5, 25, by = 5),
                                 99 + 26 + seq(5, 25, by = 5), 99 + 2 * 26 + seq(5, 25, by = 5))
                    axis(side = 2, at = -yaxis_at,
                         labels = c(rep(1984:2016, 3), rep(1991:2016, 3))[yaxis_at],
                         las = 1)
                    mtext(text = rep(c("other", "7-18", "0-6"), 2), side = 2,
                          at = -c(17, 17 + 33, 17 + 66, 99 + 13.6, 
                                  99 + 13.6 + 26, 99 + 13.6 + 2 * 26),
                          line = 4, las = 1, cex = 0.75)
                    mtext(text = c("West", "East"), side = 2,
                          at = -c(17 + 33, 99 + 13.6 + 26),
                          line = 7.3, # las = 1, 
                          cex = 0.75)
                    sapply(seq_along(hlines), 
                           function(h) segments(x0 = 0, y0 = hlines[h], x1 = hxpos[h], 
                                                y1 = hlines[h], lwd = 0.6, xpd = TRUE))
                  }
                  if (i == 1 & j %in% 1:2) {
                    segments(x0 = 0, y0 = -180.5, x1 = -1.1, 
                             y1 = -180.5, xpd = TRUE)
                  }
                }
         ))
at_y <- 177 / 2 * 1.1 # 3 + 3 * (N_main_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{f}_{West\\_East, c\\_age, year}$"), at = at_y, line = 32.5, side = 2)
mtext(text = TeX("$empCR(\\hat{f}_{West\\_East, c\\_age, year})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# # par(mfrow = c(4, 3))
# sapply(seq_along(empirical_coverage_predictions_Vc), 
#        function(i) sapply(seq_along(empirical_coverage_predictions_Vc[[i]]), 
#                           function(j) {
#                             plot(1:177, empirical_coverage_predictions_Vc[[i]][[j]], 
#                                  main = paste0("G = ", G[i], ", N = ", N[j]), ylim = c(0, 1))
#                             abline(h = 0.95, col = "green")
#                             }))

### Partial Effects (simultaneous)
coverage_effects_smt <- list()
empirical_coverage_effects_smt_Vc <- list()
empirical_coverage_effects_smt_Vp <- list()
for (i in seq_along(G)) {
  coverage_effects_smt[[i]] <- list()
  empirical_coverage_effects_smt_Vc[[i]] <- list()
  empirical_coverage_effects_smt_Vp[[i]] <- list()
  for (j in seq_along(N)) {
    coverage_effects_smt[[i]][[j]] <- list.files(path = path, 
                                                 pattern = paste0("coverage_effects_smt.*N_",
                                                                  N[j], "_G_", G[i], ".rds"), 
                                                 full.names = TRUE)
    coverage_effects_smt[[i]][[j]] <- lapply(coverage_effects_smt[[i]][[j]], readRDS)
    empirical_coverage_effects_smt_Vc[[i]][[j]] <- lapply(coverage_effects_smt[[i]][[j]], function(coverage) coverage[, 1]) # 1st column: coverage_Vc
    empirical_coverage_effects_smt_Vc[[i]][[j]] <- apply(matrix(unlist(empirical_coverage_effects_smt_Vc[[i]][[j]]), 
                                                           byrow = TRUE, 
                                                           nrow = length(empirical_coverage_effects_smt_Vc[[i]][[j]])), 
                                                    2, mean)
    empirical_coverage_effects_smt_Vp[[i]][[j]] <- lapply(coverage_effects_smt[[i]][[j]], function(coverage) coverage[, 2]) # 2nd column: coverage_Vp
    empirical_coverage_effects_smt_Vp[[i]][[j]] <- apply(matrix(unlist(empirical_coverage_effects_smt_Vp[[i]][[j]]), 
                                                           byrow = TRUE, 
                                                           nrow = length(empirical_coverage_effects_smt_Vp[[i]][[j]])), 
                                                    2, mean)
  }
}

# # Vc
# sapply(empirical_coverage_effects_smt_Vc, 
#        function(coverage) {
#          plot(seq_along(coverage), coverage, ylim = c(0, 1))
#          abline(h = 0.95, col = "green")
#        })
# 
# 
# # par(mfrow = c(4, 3))
# sapply(seq_along(empirical_coverage_predictions_Vp), 
#        function(i) sapply(seq_along(empirical_coverage_predictions_Vp[[i]]), 
#                           function(j) {
#                             plot(1:177, empirical_coverage_predictions_Vp[[i]][[j]], 
#                                  main = paste0("G = ", G[i], ", N = ", N[j]), ylim = c(0, 1))
#                             abline(h = 0.95, col = "green")
#                           }))

# params_main <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_main_effects,
#                                           byrow = FALSE, up = 3, le_ri = 8.5)
pdf("./Images/coverage_Vc_smt_main_with_pred.pdf", height = 4, width = 6.5)
layout(matrix(1:(length(G) * N_main_effects), ncol = length(G)),
       heights = c(1, rep(0.47, N_main_effects - 2), 1),
       widths = c(1, rep(0.56, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_main_effects),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_main$mar[[param_ind]])
                  if (k == 1) {
                  boxplot(empirical_coverage_predictions_Vc[[i]][length(N):1],
                          # main = main,
                          horizontal = TRUE, lwd = 0.6,
                          ylim = c(0, 1),
                          xaxt = params_main$xaxt[param_ind], yaxt = "n")
                  } else {
                    plot(sapply(empirical_coverage_effects_smt_Vc[[i]],
                                function(coverage) coverage[main_effects[k] - 1]),
                         seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                         xaxt = "n", # params_main$xaxt[param_ind],
                         yaxt =  "n")
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N),
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[main_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_main$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3 + 3 * (N_main_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

main_effects_wo_pred <- c(1:5)[-4]
N_main_effects_wo_pred <- length(main_effects_wo_pred)
params_main_wo_pred <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_main_effects_wo_pred,
                                          byrow = FALSE, up = 3, le_ri = 8.5)
pdf("./Images/coverage_Vc_smt_main.pdf", height = 3.5, width = 6.5)
layout(matrix(1:(length(G) * N_main_effects_wo_pred), ncol = length(G)),
       heights = c(1, rep(0.47, N_main_effects_wo_pred - 2), 1),
       widths = c(1, rep(0.56, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_main_effects_wo_pred),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects_wo_pred + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_main_wo_pred$mar[[param_ind]])
                    plot(sapply(empirical_coverage_effects_smt_Vc[[i]],
                                function(coverage) coverage[main_effects_wo_pred[k]]),
                         seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                         xaxt = "n", # params_main$xaxt[param_ind],
                         yaxt =  "n")
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N),
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[main_effects_wo_pred[k] + 1], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_main_wo_pred$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3 + 3 * (N_main_effects_wo_pred / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# Vp
pdf("./Images/coverage_Vp_smt_main_with_pred.pdf", height = 4, width = 6.5)
layout(matrix(1:(length(G) * N_main_effects), ncol = length(G)),
       heights = c(1, rep(0.47, N_main_effects - 2), 1),
       widths = c(1, rep(0.56, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_main_effects),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_main$mar[[param_ind]])
                  if (k == 1) {
                    boxplot(empirical_coverage_predictions_Vp[[i]][length(N):1],
                            horizontal = TRUE, lwd = 0.6,
                            ylim = c(0, 1),
                            xaxt = params_main$xaxt[param_ind], yaxt = "n")
                  } else {
                    plot(sapply(empirical_coverage_effects_smt_Vp[[i]],
                                function(coverage) coverage[main_effects[k] - 1]),
                         seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                         xaxt = "n", # params_main$xaxt[param_ind], 
                         yaxt =  "n")
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[main_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_main$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3 + 3 * (N_main_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# main_effects_wo_pred <- c(1:5)[-4]
# N_main_effects_wo_pred <- length(main_effects_wo_pred)
# params_main_wo_pred <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_main_effects_wo_pred,
#                                                  byrow = FALSE, up = 3, le_ri = 8.5)
pdf("./Images/coverage_Vp_smt_main.pdf", height = 3.5, width = 6.5)
layout(matrix(1:(length(G) * N_main_effects_wo_pred), ncol = length(G)),
       heights = c(1, rep(0.47, N_main_effects_wo_pred - 2), 1),
       widths = c(1, rep(0.56, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_main_effects_wo_pred),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects_wo_pred + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_main_wo_pred$mar[[param_ind]])
                  plot(sapply(empirical_coverage_effects_smt_Vp[[i]],
                              function(coverage) coverage[main_effects_wo_pred[k]]),
                       seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                       xaxt = "n", # params_main$xaxt[param_ind],
                       yaxt =  "n")
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N),
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[main_effects_wo_pred[k] + 1], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_main_wo_pred$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3 + 3 * (N_main_effects_wo_pred / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

### all partial effects
# params_all_MSE <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_all_effects,
#                                              byrow = FALSE, up = 3, le_ri = 13.5)
pdf("./Images/coverage_Vc_smt_all_with_pred.pdf", height = 6.3, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_all_effects),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_all$mar[[param_ind]])
                  if (k == 1) {
                    boxplot(empirical_coverage_predictions_Vc[[i]][length(N):1],
                            # main = main,
                            horizontal = TRUE, lwd = 0.6,
                            ylim = c(0, 1),
                            xaxt = params_all$xaxt[param_ind], yaxt = "n")
                  } else {
                    plot(sapply(empirical_coverage_effects_smt_Vc[[i]],
                                function(coverage) coverage[all_effects[k] - 1]),
                         seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                         xaxt = "n", # params_all$xaxt[param_ind], 
                         yaxt =  "n")
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[all_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3 + 3 * (N_all_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()


all_effects_wo_pred <- all_effects[-1] - 1
N_all_effects_wo_pred <- length(all_effects_wo_pred)
params_all_wo_pred <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_all_effects_wo_pred,
                                                  byrow = FALSE, up = 3, le_ri = 13.5)
pdf("./Images/coverage_Vc_smt_all.pdf", height = 5.8, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects_wo_pred), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects_wo_pred - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_all_effects_wo_pred),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects_wo_pred + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_all_wo_pred$mar[[param_ind]])
                  plot(sapply(empirical_coverage_effects_smt_Vc[[i]],
                              function(coverage) coverage[all_effects_wo_pred[k]]),
                       seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                       xaxt = "n", # params_all$xaxt[param_ind],
                       yaxt =  "n")
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N),
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[all_effects_wo_pred[k] + 1], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all_wo_pred$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3.5 + 3 * (N_all_effects_wo_pred / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# Vp
pdf("./Images/coverage_Vp_smt_all_with_pred.pdf", height = 6.3, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_all_effects),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_all$mar[[param_ind]])
                  if (k == 1) {
                    boxplot(empirical_coverage_predictions_Vp[[i]][length(N):1],
                            horizontal = TRUE, lwd = 0.6,
                            ylim = c(0, 1),
                            xaxt = params_all$xaxt[param_ind], yaxt = "n")
                  } else {
                    plot(sapply(empirical_coverage_effects_smt_Vp[[i]],
                                function(coverage) coverage[all_effects[k] - 1]),
                         seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                         xaxt = "n", # params_all$xaxt[param_ind], 
                         yaxt =  "n")
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[all_effects[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3 + 3 * (N_all_effects / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# # without predictions
# all_effects_wo_pred <- all_effects[-1] - 1
# N_all_effects_wo_pred <- length(all_effects_wo_pred)
# params_all_wo_pred <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_all_effects_wo_pred,
#                                                  byrow = FALSE, up = 3, le_ri = 13.5)
pdf("./Images/coverage_Vp_smt_all.pdf", height = 5.8, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects_wo_pred), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects_wo_pred - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_len(N_all_effects_wo_pred),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects_wo_pred + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_all_wo_pred$mar[[param_ind]])
                  plot(sapply(empirical_coverage_effects_smt_Vp[[i]],
                              function(coverage) coverage[all_effects_wo_pred[k]]),
                       seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                       xaxt = "n", # params_all$xaxt[param_ind],
                       yaxt =  "n")
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N),
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names[all_effects_wo_pred[k] + 1], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all_wo_pred$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3.5 + 3 * (N_all_effects_wo_pred / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1 * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# coverage_predictions_Vp <- list()
# empirical_coverage_predictions_Vp <- list()
# for (i in seq_along(N)) {
#   coverage_predictions_Vp[[i]] <- list.files(path = path, pattern = paste0("coverage_predictions_Vp_.*N_", N[i], "_G_", G, ".rds"), full.names = TRUE)
#   coverage_predictions_Vp[[i]] <- lapply(coverage_predictions_Vp[[i]], readRDS)
#   empirical_coverage_predictions_Vp[[i]] <- apply(matrix(unlist(coverage_predictions_Vp[[i]]), 
#                                                          byrow = TRUE, 
#                                                          nrow = length(coverage_predictions_Vp[[i]])), 
#                                                   2, mean)
# }
# sapply(empirical_coverage_predictions_Vp, 
#        function(coverage) {
#          plot(1:177, coverage, ylim = c(0, 1))
#          abline(h = 0.95, col = "green")
#        })

### Partial Effects (pointwise)
# Vc
coverage_effects_pw_Vc <- list()
empirical_coverage_effects_pw_Vc <- list()
for (i in seq_along(G)) {
  coverage_effects_pw_Vc[[i]] <- list()
  empirical_coverage_effects_pw_Vc[[i]] <- list()
  for (j in seq_along(N)) {
    coverage_effects_pw_Vc[[i]][[j]] <- list.files(path = path,
                                                    pattern = paste0("coverage_effects_pw_Vc_.*N_",
                                                                     N[j], "_G_", G[i], ".rds"),
                                                    full.names = TRUE)
    coverage_effects_pw_Vc[[i]][[j]] <- lapply(coverage_effects_pw_Vc[[i]][[j]], readRDS)
    cov <- coverage_effects_pw_Vc[[i]][[j]]
    empirical_coverage_effects_pw_Vc[[i]][[j]] <- lapply(seq_along(cov[[1]]),
                                                         function(n)  {
                                                           combined <- sapply(cov, function(c) c[[n]])
                                                           if (is.vector(combined)) {
                                                             mean(combined)
                                                           } else {
                                                             apply(combined, 1, mean, na.rm = TRUE)
                                                           }
                                                         })
    # empirical_coverage_effects_pw_Vc[[i]][[j]] <- apply(matrix(unlist(coverage_effects_pw_Vc[[i]][[j]]),
    #                                                             byrow = TRUE,
    #                                                             nrow = length(coverage_effects_pw_Vc[[i]][[j]])),
    #                                                      2, mean)
  }
}

# # Sollten simultan und punktweise für Intercept nicht das gleiche ergeben?
# # Normalerweise ja, aber wir haben punktweise umkodiert, sodass 1991 Referenz
# # ist, simultan nicht
# test <- lapply(1:200, function(i) readRDS(paste0(path, "coverage_effects_pw_Vp_", i, "_N_",
#                                                  N[1], "_G_", G[1], ".rds")))
# mean(sapply(test, function(t) t[[1]]))
# statistic <- lapply(1:200, function(i) readRDS(paste0(path, "statistic_effects_pw_Vp_", i, "_N_",
#                                                       N[1], "_G_", G[1], ".rds")))

library(latex2exp)
effect_names_pw <- c(TeX("$\\hat{\\f}$"), TeX("$\\hat{\\beta}_0$"), 
                     TeX("$\\hat{\\beta}_{East}$"), TeX("$\\hat{\\beta}_{7-18}$"), 
                     TeX("$\\hat{\\beta}_{0-6}$"), TeX("$\\hat{\\beta}_{7-18, East}$"),
                     TeX("$\\hat{\\beta}_{0-6, East}$"),
                     TeX("$\\hat{g}(year)$"), TeX("$\\hat{g}_{East}(year)$"),
                     TeX("$\\hat{g}_{7-18}(year)$"), TeX("$\\hat{g}_{0-6}(year)$"), 
                     TeX("$\\hat{g}_{7-18, East}(year)$"),TeX("$\\hat{g}_{0-6, East}(year)$"))

# main effects (main part of paper)
# Adjust order of c_age effects (from 7-18, 0-6 to 0-6, 7-18)
main_effects_pw <- c(1:3, 5:4, 8)
N_main_effects_pw <- length(main_effects_pw)

params_main_pw <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_main_effects_pw,
                                             byrow = FALSE, up = 3, le_ri = 8.5)
pdf("./Images/coverage_Vc_pw_main.pdf", height = 4.5, width = 6.5)
layout(matrix(1:(length(G) * N_main_effects_pw), ncol = length(G)),
       heights = c(1, rep(0.47, N_main_effects_pw - 2), 1),
       widths = c(1, rep(0.56, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_along(main_effects_pw),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects_pw + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_main_pw$mar[[param_ind]])
                  if (k == 1) {
                    boxplot(empirical_coverage_predictions_Vc[[i]][length(N):1],
                            # main = main,
                            horizontal = TRUE, lwd = 0.6,
                            ylim = c(0, 1),
                            xaxt = "n", yaxt = "n")
                  } else {
                    cov_effects <- sapply(empirical_coverage_effects_pw_Vc[[i]],
                                          function(coverage) coverage[main_effects_pw[k] - 1])
                    if (length(cov_effects[[1]]) == 1) {
                      plot(unlist(cov_effects),
                           seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                           xaxt = "n", yaxt =  "n")
                    } else {
                      boxplot(cov_effects[length(N):1], horizontal = TRUE, lwd = 0.6,
                              ylim = c(0, 1), xaxt = "n", yaxt = "n")
                    }
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names_pw[main_effects_pw[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_main_pw$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3.5 + 3 * (N_main_effects_pw / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# all effects (appendix)
# Adjust order of c_age effects (from 7-18, 0-6 to 0-6, 7-18)
all_effects_pw <- c(1:3, 5:4, 7:6, 8:9, 11:10, 13:12) # 1:13
N_all_effects_pw <- length(all_effects_pw)

params_all_pw <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_all_effects_pw,
                                            byrow = FALSE, up = 3, le_ri = 13.5)
pdf("./Images/coverage_Vc_pw_all.pdf", height = 8.3, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects_pw), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects_pw - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_along(all_effects_pw),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects_pw + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_all_pw$mar[[param_ind]])
                  if (k == 1) {
                    boxplot(empirical_coverage_predictions_Vc[[i]][length(N):1],
                            # main = main,
                            horizontal = TRUE, lwd = 0.6,
                            ylim = c(0, 1), xaxt = "n", yaxt = "n")
                  } else {
                    cov_effects <- sapply(empirical_coverage_effects_pw_Vc[[i]],
                                          function(coverage) coverage[all_effects_pw[k] - 1])
                    if (length(cov_effects[[1]]) == 1) {
                      plot(unlist(cov_effects),
                           seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                           xaxt = "n", yaxt =  "n")
                    } else {
                      boxplot(cov_effects[length(N):1], horizontal = TRUE, lwd = 0.6,
                              ylim = c(0, 1), xaxt = "n", yaxt = "n")
                    }
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names_pw[all_effects_pw[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all_pw$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3.5 + 3 * (N_all_effects_pw / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# Vp
coverage_effects_pw_Vp <- list()
empirical_coverage_effects_pw_Vp <- list()
for (i in seq_along(G)) {
  coverage_effects_pw_Vp[[i]] <- list()
  empirical_coverage_effects_pw_Vp[[i]] <- list()
  for (j in seq_along(N)) {
    coverage_effects_pw_Vp[[i]][[j]] <- list.files(path = path,
                                                   pattern = paste0("coverage_effects_pw_Vp_.*N_",
                                                                    N[j], "_G_", G[i], ".rds"),
                                                   full.names = TRUE)
    coverage_effects_pw_Vp[[i]][[j]] <- lapply(coverage_effects_pw_Vp[[i]][[j]], readRDS)
    cov <- coverage_effects_pw_Vp[[i]][[j]]
    empirical_coverage_effects_pw_Vp[[i]][[j]] <- lapply(seq_along(cov[[1]]),
                                                         function(n)  {
                                                           combined <- sapply(cov, function(c) c[[n]])
                                                           if (is.vector(combined)) {
                                                             mean(combined)
                                                           } else {
                                                             apply(combined, 1, mean, na.rm = TRUE)
                                                           }
                                                         })
    # empirical_coverage_effects_pw_Vp[[i]][[j]] <- apply(matrix(unlist(coverage_effects_pw_Vp[[i]][[j]]),
    #                                                            byrow = TRUE,
    #                                                            nrow = length(coverage_effects_pw_Vp[[i]][[j]])),
    #                                                     2, mean)
  }
}


# library(latex2exp)
# effect_names_pw <- c(TeX("$\\hat{\\f}$"), TeX("$\\hat{\\beta}_0$"), 
#                      TeX("$\\hat{\\beta}_{East}$"), TeX("$\\hat{\\beta}_{7-18}$"), 
#                      TeX("$\\hat{\\beta}_{0-6}$"), TeX("$\\hat{\\beta}_{7-18, East}$"),
#                      TeX("$\\hat{\\beta}_{0-6, East}$"),
#                      TeX("$\\hat{g}(year)$"), TeX("$\\hat{g}_{East}(year)$"),
#                      TeX("$\\hat{g}_{7-18}(year)$"), TeX("$\\hat{g}_{0-6}(year)$"), 
#                      TeX("$\\hat{g}_{7-18, East}(year)$"),TeX("$\\hat{g}_{0-6, East}(year)$"))

# main effects (main part of paper)
# # Adjust order of c_age effects (from 7-18, 0-6 to 0-6, 7-18)
# main_effects_pw <- c(1:3, 5:4, 8)
# N_main_effects_pw <- length(main_effects_pw)

# params_main_pw <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_main_effects_pw,
#                                              byrow = FALSE, up = 3, le_ri = 8.5)
pdf("./Images/coverage_Vp_pw_main.pdf", height = 4.5, width = 6.5)
layout(matrix(1:(length(G) * N_main_effects_pw), ncol = length(G)),
       heights = c(1, rep(0.47, N_main_effects_pw - 2), 1),
       widths = c(1, rep(0.56, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_along(main_effects_pw),
                function(k) {
                  param_ind <- (i - 1) * N_main_effects_pw + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_main_pw$mar[[param_ind]])
                  if (k == 1) {
                    boxplot(empirical_coverage_predictions_Vp[[i]][length(N):1],
                            # main = main,
                            horizontal = TRUE, lwd = 0.6,
                            ylim = c(0, 1),
                            xaxt = "n", yaxt = "n")
                  } else {
                    cov_effects <- sapply(empirical_coverage_effects_pw_Vp[[i]],
                                          function(coverage) coverage[main_effects_pw[k] - 1])
                    if (length(cov_effects[[1]]) == 1) {
                      plot(unlist(cov_effects),
                           seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                           xaxt = "n", yaxt =  "n")
                    } else {
                      boxplot(cov_effects[length(N):1], horizontal = TRUE, lwd = 0.6,
                              ylim = c(0, 1), xaxt = "n", yaxt = "n")
                    }
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names_pw[main_effects_pw[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_main_pw$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3.5 + 3 * (N_main_effects_pw / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# # all effects (appendix)
# # Adjust order of c_age effects (from 7-18, 0-6 to 0-6, 7-18)
# all_effects_pw <- c(1:3, 5:4, 7:6, 8:9, 11:10, 13:12) # 1:13
# N_all_effects_pw <- length(all_effects_pw)

# params_all_pw <- get_params_for_matrix_plot(n_cols = length(G), n_rows = N_all_effects_pw,
#                                              byrow = FALSE, up = 3, le_ri = 13.5)
pdf("./Images/coverage_Vp_pw_all.pdf", height = 8.3, width = 7.8)
layout(matrix(1:(length(G) * N_all_effects_pw), ncol = length(G)),
       heights = c(1, rep(0.47, N_all_effects_pw - 2), 1),
       widths = c(1, rep(0.44, length(G) - 2), 1))
sapply(seq_along(G),
       function(i)
         sapply(seq_along(all_effects_pw),
                function(k) {
                  param_ind <- (i - 1) * N_all_effects_pw + k
                  if (k == 1) {
                    main <- paste0("G = ", G[i])
                  } else {
                    main <- ""
                  }
                  par(mar = params_all_pw$mar[[param_ind]])
                  if (k == 1) {
                    boxplot(empirical_coverage_predictions_Vp[[i]][length(N):1],
                            # main = main,
                            horizontal = TRUE, lwd = 0.6,
                            ylim = c(0, 1), xaxt = "n", yaxt = "n")
                  } else {
                    cov_effects <- sapply(empirical_coverage_effects_pw_Vp[[i]],
                                          function(coverage) coverage[all_effects_pw[k] - 1])
                    if (length(cov_effects[[1]]) == 1) {
                      plot(unlist(cov_effects),
                           seq_along(N)[length(N):1], xlim = c(0, 1), xlab = "", ylab = "", ylim = c(0.5, 3.5),
                           xaxt = "n", yaxt =  "n")
                    } else {
                      boxplot(cov_effects[length(N):1], horizontal = TRUE, lwd = 0.6,
                              ylim = c(0, 1), xaxt = "n", yaxt = "n")
                    }
                  }
                  abline(v = 0.95, col = "green3", lwd = 0.8)
                  mtext(text = main, side = 3, line = 1)
                  if (i == length(G)) {
                    # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
                    mtext(side = 4, at = length(N):1, text = paste0("N = ", N), 
                          las = 1, line = 1)
                  }
                  if (i == 1) {
                    mtext(text = effect_names_pw[all_effects_pw[k]], # rownames(relMSE[[i]][[1]])[k],
                          side = 2, las = 1, line = 1)
                  }
                  if (params_all_pw$xaxt[param_ind] == "s") {
                    axis(side = 1, at = seq(0, 1, by = 0.2),
                         labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
                  }
                }
         ))
at_y <- 3.5 + 3 * (N_all_effects_pw / 2 - 1) * 1.1
at_x <- - (length(G) / 2 - 1) * 1.1
mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 33.1, side = 2)
mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
dev.off()

# # Boxplot over all effects:
# pdf("./Images/coverage_Vp_pw_boxplot.pdf", height = 4, width = 6.5)
# layout(matrix(1:(length(G) * length(N)), ncol = length(G)),
#        heights = c(1, rep(0.8, length(N) - 2), 1),
#        widths = c(1, rep(0.44, length(G) - 2), 1))
# sapply(seq_along(G),
#        function(i)
#          sapply(seq_along(N),
#                 function(j) {
#                   param_ind <- (i - 1) * length(N) + j
#                   if (j == 1) {
#                     main <- paste0("G = ", G[i])
#                   } else {
#                     main <- ""
#                   }
#                   par(mar = params_pred$mar[[param_ind]])
#                   boxplot(empirical_coverage_effects_pw_Vp[[i]][[j]],
#                             horizontal = TRUE, lwd = 0.6,
#                             ylim = c(0, 1),
#                             xaxt = "n", yaxt = "n")
#                   abline(v = 0.95, col = "green3", lwd = 0.8)
#                   mtext(text = main, side = 3, line = 1)
#                   if (i == length(G)) {
#                     # axis(side = 4, at = length(N):1, labels = paste0("N = ", N), las = 1)
#                     mtext(side = 4, at = 1, text = paste0("N = ", N[j]), 
#                           las = 1, line = 1)
#                   }
#                   # if (i == 1) {
#                   #   mtext(text = effect_names[all_effects[k]], # rownames(relMSE[[i]][[1]])[k],
#                   #         side = 2, las = 1, line = 1)
#                   # }
#                   if (params_pred$xaxt[param_ind] == "s") {
#                     axis(side = 1, at = seq(0, 1, by = 0.2),
#                          labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
#                   }
#                 }
#          ))
# at_y <- 3 + 3 * (N_all_effects / 2 - 1) * 1.1
# at_x <- - (length(G) / 2 - 1) * 1.1
# mtext(text = TeX("$\\hat{e}$"), at = at_y, line = 28.5, side = 2)
# mtext(text = TeX("$empCR(\\hat{e})$"), at = at_x, line = 3.5, side = 1)
# dev.off()

########################### Evaluate partial effects ###########################
### relMSE

# Gespeicherter relMSE betrachtet einzelne Kategorien für diskrete Kovariablen
# separat:
# relMSE_effects <- list()
# for (i in seq_along(N)) {
#   relMSE_effects[[i]] <- list.files(path = path, pattern = paste0("relMSE_effects.*N_", N[i], "_G_", G, ".rds"), full.names = TRUE)
#   relMSE_effects[[i]] <- lapply(relMSE_effects[[i]], readRDS)
#   relMSE_effects[[i]] <- matrix(unlist(relMSE_effects[[i]]),
#                                 byrow = TRUE,
#                                 nrow = length(relMSE_effects[[i]]))
# }
# par(mfrow = c(1, 2))
# sapply(relMSE_effects, boxplot)
# lapply(relMSE_effects, function(relMSE) apply(relMSE, 2, fivenum))

# Analog zum 1. Paper für diskrete Kovariablen über Kategorien mitteln:
# c_age: 3, 4
# West_East_c_age: 5, 6
# year_c_age: 9, 10
# year_West_East_c_age: 11, 12
average_categories <- function(x) {
  c(x[1:2], mean(x[3:4]), mean(x[5:6]), x[7:8], mean(x[9:10]), mean(x[11:12]))
}

denominator_effects <- readRDS(paste0(path, "denominator_effects_relMSE.rds"))
denominator_effects_agg <- average_categories(denominator_effects)

MSE_effects <- list()
relMSE_effects <- list()
for (i in seq_along(N)) {
  MSE_effects[[i]] <- list.files(path = path, pattern = paste0("^MSE_effects.*N_", N[i], "_G_", G, ".rds"), full.names = TRUE)
  MSE_effects[[i]] <- lapply(MSE_effects[[i]], readRDS)
  MSE_effects[[i]] <- lapply(MSE_effects[[i]], average_categories)
  relMSE_effects[[i]] <- lapply(MSE_effects[[i]], function(MSE) MSE / denominator_effects_agg)
  relMSE_effects[[i]] <- matrix(unlist(relMSE_effects[[i]]), byrow = TRUE, 
                             nrow = length(relMSE_effects[[i]]))
}

par(mfrow = c(1, 2))
sapply(relMSE_effects, boxplot)
lapply(relMSE_effects, function(relMSE) apply(relMSE, 2, fivenum))

### Simultaneous Coverage
coverage_effects_smt <- list()
empirical_coverage_effects_smt_Vc <- list()
empirical_coverage_effects_smt_Vp <- list()
for (i in seq_along(N)) {
  coverage_effects_smt[[i]] <- list.files(path = path, pattern = paste0("coverage_effects_smt.*N_", N[i], "_G_", G, ".rds"), full.names = TRUE)
  coverage_effects_smt[[i]] <- lapply(coverage_effects_smt[[i]], readRDS)
  empirical_coverage_effects_smt_Vc[[i]] <- lapply(coverage_effects_smt[[i]], function(coverage) coverage[, 1]) # 1st column: coverage_Vc
  empirical_coverage_effects_smt_Vc[[i]] <- apply(matrix(unlist(empirical_coverage_effects_smt_Vc[[i]]), 
                                                         byrow = TRUE, 
                                                         nrow = length(empirical_coverage_effects_smt_Vc[[i]])), 
                                                  2, mean)
  empirical_coverage_effects_smt_Vp[[i]] <- lapply(coverage_effects_smt[[i]], function(coverage) coverage[, 2]) # 2nd column: coverage_Vp
  empirical_coverage_effects_smt_Vp[[i]] <- apply(matrix(unlist(empirical_coverage_effects_smt_Vp[[i]]), 
                                                         byrow = TRUE, 
                                                         nrow = length(empirical_coverage_effects_smt_Vp[[i]])), 
                                                  2, mean)
}
# Vc
sapply(empirical_coverage_effects_smt_Vc, 
       function(coverage) {
         plot(seq_along(coverage), coverage, ylim = c(0, 1))
         abline(h = 0.95, col = "green")
       })

# Vp
sapply(empirical_coverage_effects_smt_Vp, 
       function(coverage) {
         plot(seq_along(coverage), coverage, ylim = c(0, 1))
         abline(h = 0.95, col = "green")
       })

### Pointwise Coverage (Vc)
coverage_effects_pw_Vc <- list()
empirical_coverage_effects_pw_Vc <- list()
for (i in seq_along(N)) {
  coverage_effects_pw_Vc[[i]] <- list.files(path = path, pattern = paste0("coverage_effects_pw_Vc_.*N_", N[i], "_G_", G, ".rds"), full.names = TRUE)
  coverage_effects_pw_Vc[[i]] <- lapply(coverage_effects_pw_Vc[[i]], readRDS)
  coverage_effects_pw_Vc[[i]] <- matrix(unlist(coverage_effects_pw_Vc[[i]]), byrow = TRUE, 
                                    nrow = length(coverage_effects_pw_Vc[[i]]))
  empirical_coverage_effects_pw_Vc[[i]] <- apply(coverage_effects_pw_Vc[[i]], 2, mean)
}
sapply(empirical_coverage_effects_pw_Vc, 
       function(coverage) {
         plot(seq_along(coverage), coverage, ylim = c(0, 1))
         abline(h = 0.95, col = "green")
       })

indices <- list(intercepts = 1:6, year = 6 + 1:33, year_east = 6 + 33 + 1:26,
                year_1 = 6 + 33 + 26 + 1:33, year_2 = 6 + 2 * 33 + 26 + 1:33,
                year_east_1 = 6 + 3 * 33 + 26 + 1:26, year_east_2 = 6 + 3 * 33 + 2 * 26 + 1:26)
names <- list(c(intercepts = "Intercept", "East", "1", "2", "East1", "East2"),
              year = 1984:2016, year_east = 1991:2016, year_1 = 1984:2016,
              year_2 = 1984:2016, year_east_1 = 1991:2016, year_east_2 = 1991:2016)


sapply(seq_along(indices), 
       function(i) {
         sapply(empirical_coverage_effects_pw_Vc, 
                function(coverage) {
                  plot(indices[[i]], coverage[indices[[i]]], ylim = c(0, 1), xlab = "year", 
                       xaxt = "n", main = names(indices)[i])
                  axis(1, at = indices[[i]], labels = names[[i]])
                  abline(h = 0.95, col = "green")
                })
       })


### Pointwise Coverage (Vp)
coverage_effects_pw_Vp <- list()
empirical_coverage_effects_pw_Vp <- list()
for (i in seq_along(N)) {
  coverage_effects_pw_Vp[[i]] <- list.files(path = path, pattern = paste0("coverage_effects_pw_Vp_.*N_", N[i], "_G_", G, ".rds"), full.names = TRUE)
  coverage_effects_pw_Vp[[i]] <- lapply(coverage_effects_pw_Vp[[i]], readRDS)
  coverage_effects_pw_Vp[[i]] <- matrix(unlist(coverage_effects_pw_Vp[[i]]), byrow = TRUE, 
                                        nrow = length(coverage_effects_pw_Vp[[i]]))
  empirical_coverage_effects_pw_Vp[[i]] <- apply(coverage_effects_pw_Vp[[i]], 2, mean)
}
sapply(empirical_coverage_effects_pw_Vp, 
       function(coverage) {
         plot(seq_along(coverage), coverage, ylim = c(0, 1))
         abline(h = 0.95, col = "green")
       })

indices <- list(intercepts = 1:6, year = 6 + 1:33, year_east = 6 + 33 + 1:26,
                year_1 = 6 + 33 + 26 + 1:33, year_2 = 6 + 2 * 33 + 26 + 1:33,
                year_east_1 = 6 + 3 * 33 + 26 + 1:26, year_east_2 = 6 + 3 * 33 + 2 * 26 + 1:26)
names <- list(c(intercepts = "Intercept", "East", "1", "2", "East1", "East2"),
              year = 1984:2016, year_east = 1991:2016, year_1 = 1984:2016,
              year_2 = 1984:2016, year_east_1 = 1991:2016, year_east_2 = 1991:2016)


sapply(seq_along(indices), 
       function(i) {
         sapply(empirical_coverage_effects_pw_Vp, 
                function(coverage) {
                  plot(indices[[i]], coverage[indices[[i]]], ylim = c(0, 1), xlab = "year", 
                       xaxt = "n", main = names(indices)[i])
                  axis(1, at = indices[[i]], labels = names[[i]])
                  abline(h = 0.95, col = "green")
                })
       })

