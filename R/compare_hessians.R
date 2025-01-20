# this script is meant to validate if the gam package uses the correct
# covariance matrix or if the covariance matrix has to be corrected
library(matlib)

cov_mat_est <- model_soep$Ve
rV <- model_soep$rV

predicted_values <- model_soep$linear.predictors
lp_matrix <- predict(model_soep, type = "lpmatrix")
smooth <- model_soep$smooth

delta <- dta_est$Delta
offset <- model_soep$offset
# now  reconstruct the penalty matrix

n_smooth <- n_smooths(model_soep)
n_params <- length(model_soep$coefficients)

n_zero_cols_start <- 0
for (h in 1:n_smooth) {
   penalty_block <- smooth[[h]]$S[[1]]
   n_cols_penalty_block <- ncol(penalty_block)
   n_row_penalty_block <- nrow(penalty_block)
   print(paste0(n_row_penalty_block, " x ", n_cols_penalty_block ))
   n_zero_cols_end <- n_params - n_zero_cols_start - n_cols_penalty_block
   zero_mat_end <- matrix(0,
                          n_row_penalty_block,
                          n_zero_cols_end)
   if (h == 1) {
     penalty_mat <- cbind(penalty_block, zero_mat_end)
   }
   else {
     zero_mat_start <- matrix(0,
                              n_row_penalty_block,
                              n_zero_cols_start)
     penalty_mat_rows <- cbind(zero_mat_start, penalty_block, zero_mat_end)
     penalty_mat <- rbind(penalty_mat, penalty_mat_rows)
   }

   n_zero_cols_start <- n_zero_cols_start + n_cols_penalty_block
}

penalty_mat_intercepts <- matrix(0, n_params - nrow(penalty_mat), n_params)
penalty_mat <- rbind(penalty_mat, penalty_mat_intercepts)

## now I can calculate the estimated hessian matrix
hessian_estimated <- MASS::ginv(cov_mat_est) - penalty_mat
inv(cov_mat_est)

## now I calculate the hessian which is expected when using a poisson model "naively"

## now I calculate the hessian of a multinomial model

## comparison of hessians and corresponding covariance matrices

