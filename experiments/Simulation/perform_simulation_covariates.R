grid <- seq(0, 1, length.out = 100)
exp_sinus <- function(x) {
  exp(sin(x  *  3.14))
}
exp_sinus_values <- exp_sinus(grid)
integrated_exp_sinus <- integrate(exp_sinus, 0, 1)
sinus_density <- exp_sinus(grid) / integrate(exp_sinus, 0, 1)$value
plot(grid, sinus_density, type = "l")
plot(exp(sin(grid * 3.14)))
beta_distribution <- dbeta(grid, 0.1, 3)
plot(grid, beta_distribution, type = "l")
plot(grid, f_0_line(grid, 0, 1), type = "l")

