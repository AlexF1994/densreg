# get_params_for_matrix_plot computes mar and x-/yaxt values for arranging plots
# in a matrix. Arguments:
# - n_cols, n_rows: number of columns/rows of the matrix
# - up_down: mar-value used for top/bottom space (first/last row of plot matrix)
# - le_ri: mar-value used for left/right space (first/last column of plot matrix)
# - byrow: 
get_params_for_matrix_plot <- function(n_cols = 3, n_rows = 4, up_down = 4.7, 
                                       le_ri = 4.6, byrow = TRUE, ri = 1, up = 1) {
  if (byrow) {
    if (n_cols >= 2) {
      mar_row <- c(list(c(0, le_ri, 0, 0)), rep(list(c(0, 0, 0, 0)), n_cols - 2), list(c(0, 0, 0, le_ri)))
    } else {
      mar_row <- list(c(0, le_ri, 0, ri)) # 1 column
    }
    if (n_rows >= 3) {
      mar_inner <- rep(mar_row, n_rows - 2)
    } else {
      mar_inner <- NULL
    }
    if (n_rows >= 2) {
      if (n_cols >= 2) {
        mar <- c(list(c(0, le_ri, up_down, 0)), rep(list(c(0, 0, up_down, 0)), n_cols - 2), list(c(0, 0, up_down, le_ri)),
                 mar_inner,
                 list(c(up_down, le_ri, 0, 0)), rep(list(c(up_down, 0, 0, 0)), n_cols - 2), list(c(up_down, 0, 0, le_ri)))
      } else {
        mar <- c(list(c(0, le_ri, up_down, ri)), mar_inner, list(c(up_down, le_ri, 0, ri)))
      }
    } else {
      if (n_cols >= 2) {
        mar <- c(list(c(up_down, le_ri, 0, 0)), rep(list(c(up_down, 0, 0, 0)), n_cols - 2), list(c(up_down, 0, 0, le_ri)))
      } else {
        mar <- list(c(up_down, le_ri, up, ri))
      }
    }
    xaxt <- rep(c(rep("n", n_rows - 1),  "s"), each = n_cols) # c(rep("n", (n_rows - 1) * 3), c("n", "s", "n"))
    yaxt <- rep(c("s", rep("n", n_cols - 1)), n_rows)
  } else {
    if (n_rows >= 2) {
      mar_col <- c(list(c(0, 0, up_down, 0)), rep(list(c(0, 0, 0, 0)), n_rows - 2), list(c(up_down, 0, 0, 0)))
    } else {
      mar_col <- list(c(up_down, 0, up, 0)) # 1 row
    }
    if (n_cols >= 3) {
      mar_inner <- rep(mar_col, n_cols - 2)
    } else {
      mar_inner <- NULL
    }
    if (n_cols >= 2) {
      if (n_rows >= 2) {
        mar <- c(list(c(0, le_ri, up_down, 0)), rep(list(c(0, le_ri, 0, 0)), n_rows - 2), list(c(up_down, le_ri, 0, 0)),
                 mar_inner,
                 list(c(0, 0, up_down, le_ri)), rep(list(c(0, 0, 0, le_ri)), n_rows - 2), list(c(up_down, 0, 0, le_ri)))
      } else {
        mar <- c(list(c(up_down, le_ri, up, 0)), mar_inner, list(c(up_down, 0, up, le_ri)))
      }
    } else {
      # mar <- list(c(up_down, le_ri, 0, 0), c(up_down, 0, 0, 0), c(up_down, 0, 0, le_ri))
      if (n_cols >= 2) {
        mar <- c(list(c(up_down, le_ri, 0, 0)), rep(list(c(up_down, 0, 0, 0)), n_cols - 2), list(c(up_down, 0, 0, le_ri)))
      } else {
        mar <- list(c(up_down, le_ri, up, ri))
      }
    }
    xaxt <- rep(c(rep("n", n_rows - 1),  "s"), n_cols) # c(rep("n", (n_rows - 1) * 3), c("n", "s", "n"))
    yaxt <- rep(c("s", rep("n", n_cols - 1)), each = n_rows)
  }
  return(list(mar = mar, xaxt = xaxt, yaxt = yaxt))
}
