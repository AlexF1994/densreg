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
features <- c("age_woman")

dta_binned = bin_features(dta, features, 20)

dta_grouped <- dta_binned[, group_id := .GRP,
                          keyby = .(age_woman,
                                    child_group,
                                    syear
                          )][, .(share,
                                 age_woman,
                                 child_group,
                                 syear,
                                 group_id)]

dta_est <- construct_dataset(dta_grouped,
                             response = list("share"),
                             step_size = step_size
)

poisson_mod_ind <- gam(counts ~ s(share, bs = "md",  by = as.factor(child_group), m=c(2,2))
                       + ti(share, age_woman, bs = c("md","ps"), m = list(c(2,2), c(2,2)))
                       + ti(share, syear, bs = c("md","ps"), m = list(c(2,2), c(2,2)))
                       + as.factor(group_id)
                       - 1,
                       data = dta_est,
                       method = "REML",
                       family = poisson())

plot(poisson_mod_ind)

# general intercept
poisson_mod_gen <- gam(counts ~ s(share, bs = "md",  by = as.factor(child_group), m=c(2,2))
                       + ti(share, age_woman, bs = c("md", "ps"), m = list(c(2,2), c(2,2)))
                       + s(age_woman, m=c(2,2)),
                       data = dta_est,
                       method = "REML",
                       family = poisson())

plot(poisson_mod_gen)
