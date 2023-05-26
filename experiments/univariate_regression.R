## This code implements the poisson density regression approach for the SOEP dataset
library(data.table)
library(mgcv)
library(ggplot2)
library(mgcViz)
library(gridExtra)
library(tidyr)
library(ggpubr)
library(compositions)
library(dplyr)
library(tidymv)
source("./experiments/data_prep.R")
source("./experiments/mixed_density_smooth.R")

# read in data to data.table
dta <- as.data.table(readRDS("C:/Users/fottneal/Documents/MA_Alexander_Fottner/income_share_data.rds"))

# now I drop all entries where the income_share is 0 or 1 and NAs
dta <- dta[is.na(age_man) == FALSE]

## looking for NAs
dta[is.na(age_man), .N] # 0 NA
dta[is.na(age_woman), .N] # 0 NA
dta[is.na(child_group), .N] # 0 NA
dta[is.na(syear), .N] # 0 NA

step_size <- 0.01
n_bins <- 1/ step_size


# univariate density estimation and hist plotting
grid_hist <- seq(from = 0, to = 1, by = step_size)
hist_dens <- hist(dta[share > 0 & share < 1]$share, breaks = grid_hist, right = TRUE)
data_d <-dta[!(share > 0 & share < 1)][, .N, by=share]
colnames(data_d) <- c("share", "counts")


counts_dens <- hist_dens$counts
mids_dens <- hist_dens$mids
dta_dens_c <- as.data.table(cbind(counts_dens, mids_dens))
colnames(dta_dens_c) <- c("counts", "share")
dta_dens <- rbind(dta_dens_c, data_d)[order(share)]

poisson_mod_dens <- gam(counts ~ s(share, bs = "md", m=c(2,2)),
                        data = dta_dens,
                        method = "REML",
                        family = poisson())

dens_unnorm <- predict(poisson_mod_dens, type = "response")
gg_data_uni <- data.table(share = dta_dens$share, predictions = dens_unnorm)

# plot histogram using ggplot
ggplot(dta, aes(x = share)) +
  geom_histogram(bins = n_bins) +
  labs(x = "share of woman's income", y = "number of observations")

# plot univariate density using ggplot
ggplot(dta, aes(x = share)) +
  geom_histogram(bins = n_bins) +
  labs(x = "share of woman's income", y = "number of observations") +
  geom_line(data = gg_data_uni,
            aes(x = share, y = predictions),
            color = "red",
            size = 1)

