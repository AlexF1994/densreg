### The following code checks for significance of all estimated effects

################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/../.."))

# Package to estimate GAMs (like our Poisson model)
library(mgcv)

# load (amongst others) function to check for significance
source("./experiments/Eva_SOEP/help_functions.R")
# load self-written smoother for mixed reference measure
source("./R/mixed_density_smooth.R")

# model was estimated in estimate_model.R
model_soep <- readRDS("./experiments/Eva_SOEP/Objects/model_all_covariates_but_region.rds")
# estimated effects were computed in plot_estimated_effects.R
effects <- readRDS("./experiments/Eva_SOEP/Objects/estimated_effects.rds")

################################################################################
############################ Check for significance ############################
################################################################################

# The smooth terms contained in model_soep$smooth[[i]] do not correspond to the
# effects of the model in our paper directly, since we use a different coding.
# To check for the significance of the (coding transformed) effects, we have to
# combine different smooth terms properly (possibly evaluated at specific
# covariate values), see Appendix ??? of the paper. For better tracability, here
# is the list of effects contained in model_soep$smooth[[i]]:
# i = 1:  beta_0; smooth intercept for reference West Germany, c_age other
# i = 2:  beta_{East}; smooth effect for East Germany
# i = 3:  beta_{2}; smooth effects for c_age 2
# i = 4:  beta_{1}; smooth effects for c_age 1
# i = 5:  beta_{2, West}; smooth interaction of c_age 2 and West Germany
# i = 6:  beta_{1, West}; smooth interaction of c_age 1 and West Germany
# i = 7:  beta_{3, East}; smooth interaction of c_age 3 and East Germany
# i = 8:  beta_{2, East}; smooth interaction of c_age 2 and East Germany
# i = 9:  beta_{1, East}; smooth interaction of c_age 1 and East Germany
# i = 10: g_0(year); smooth effect for year
# i = 11: g_{East}(year); smooth interaction of year and East Germany
# i = 12: g_{2}(year); smooth interaction of year and c_age 2
# i = 13: g_{1}(year); smooth interaction of year and c_age 1
# i = 14: g_{2, West}(year); smooth interaction of year, c_age 2 and West Germany
# i = 15: g_{1, West}(year); smooth interaction of year, c_age 1 and West Germany
# i = 16: g_{3, East}(year); smooth interaction of year, c_age 3 and East Germany
# i = 17: g_{2, East}(year); smooth interaction of year, c_age 2 and East Germany
# i = 18: g_{1, East}(year); smooth interaction of year, c_age 1 and East Germany

# Note that all effects corresponding to the reference category are automatically
# 0 and do not appear here! (E.g., beta_{West}, beta_{3, West})

### intercept = beta_0 + g_0(1991)
significance_intercept <- check_significance(model_soep, n_smooth = c(1, 10),
                                             n_x = list(NULL, 8), effect = effects$intercept) # 1991 is 8-th observed year
significance_intercept$significant

### West_East = beta_{West_East} + beta{3, West_East}
###           + g_{West_East}(1991) + g_{3, West_East}(1991)
significance_East <- check_significance(model_soep, n_smooth = c(2, 7, 11, 16),
                                        n_x = list(NULL, NULL, 1, 1), # in East, 1991 is the first observed year
                                        effect = effects$West_East[2, ])
significance_East$significant

### c_age = beta_{c_age} + beta_{c_age, West}
###         + g_{c_age}(1991) + g_{c_age, West}(1991)
significance_2 <- check_significance(model_soep, n_smooth = c(3, 5, 12, 14),
                                     n_x = list(NULL, NULL, 8, 8),
                                     effect = effects$c_age[2, ])
significance_2$significant

significance_1 <- check_significance(model_soep, n_smooth = c(4, 6, 13, 15),
                                     n_x = list(NULL, NULL, 8, 8),
                                     effect = effects$c_age[1, ])
significance_1$significant

### c_age_West_East = beta_{c_age, West_East} - beta_{3, West_East} - beta_{c_age, West}
###                 + g_{c_age, West_East}(1991) - g_{3, West_East}(1991) - g_{c_age, West}(1991)
# For West, all effects are 0, for East, interaction with c_age 3 is 0
significance_2_East <- check_significance(model_soep,
                                          n_smooth = c(8, -7, -5, 17, -16, -14),
                                          n_x = list(NULL, NULL, NULL, 1, 1, 8),
                                          effect = effects$West_East_c_age[5, ])
significance_2_East$significant

significance_1_East <- check_significance(model_soep,
                                          n_smooth = c(9, -7, -6, 18, -16, -15),
                                          n_x = list(NULL, NULL, NULL, 1, 1, 8),
                                          effect = effects$West_East_c_age[4, ])
significance_1_East$significant

### year = g_0(year) - g_0(1991)
significance_year <- check_significance(model_soep, c(10, -10), n_x = list(NULL, 8),
                                        effect = effects$year_intercept)
sum(unlist(significance_year$significant), na.rm = TRUE)
which(is.na(unlist(significance_year$significant))) # reference year 1991 is NA, as expected

### year_{West_East} = g_{West_East}(year) + g_{3, West_East}(year)
###                  - g_{West_East}(1991) - g_{3, West_East}(1991)
significance_year_East <- check_significance(model_soep,
                                             n_smooth = c(11, 16, -11, -16),
                                             n_x = list(NULL, NULL, 1, 1),
                                             effect = effects$year_East[-(1:7),])
sum(unlist(significance_year_East$significant), na.rm = TRUE)
which(is.na(unlist(significance_year_East$significant))) # reference year 1991 is NA, as expected
which(unlist(significance_year_East$significant)) + 1990 # 1992-2009 are not significant, 2010-2016 are significant

### year_{c_age} = g_{c_age}(year) + g_{c_age, West}(year)
###                - g_{c_age}(1991) - g_{c_age, West}(1991)
significance_year_2 <- check_significance(model_soep,
                                          n_smooth = c(12, 14, -12, -14),
                                          n_x = list(NULL, NULL, 8, 8),
                                          effect = effects$year_2)
sum(unlist(significance_year_2$significant), na.rm = TRUE)
which(is.na(unlist(significance_year_2$significant))) # reference year 1991 is NA, as expected

significance_year_1 <- check_significance(model_soep,
                                          n_smooth = c(13, 15, -13, -15),
                                          n_x = list(NULL, NULL, 8, 8),
                                          effect = effects$year_1)
sum(unlist(significance_year_2$significant), na.rm = TRUE)
which(is.na(unlist(significance_year_2$significant))) # reference year 1991 is NA, as expected

### year_{c_age, West_East} = g_{c_age, West_East}(year) - g_{3, West_East}(year)
###                         - g_{c_age, West}(year)
###                         - g_{c_age, West_East}(1991) + g_{3, West_East}(1991)
###                         + g_{c_age, West}(1991)
# For West, all effects are 0, for East, interaction with c_age 3 is 0
significance_year_2_East <- check_significance(model_soep,
                                               n_smooth = c(17, -16, -14, -17, 16, 14),
                                               n_x = list(NULL, NULL, NULL, 1, 1, 8),
                                               effect = effects$year_East_2[-(1:7),])
sum(unlist(significance_year_2_East$significant), na.rm = TRUE)
which(is.na(unlist(significance_year_2_East$significant))) # reference year 1991 is NA, as expected
which(unlist(significance_year_2_East$significant)) + 1990 # 1992-1998 and 2011-2016 are not significant, 1999-2010 are significant

significance_year_1_East <- check_significance(model_soep,
                                               n_smooth = c(18, -16, -15, -18, 16, 15),
                                               n_x = list(NULL, NULL, NULL, 1, 1, 8),
                                               effect = effects$year_East_1[-(1:7),])
sum(unlist(significance_year_1_East$significant), na.rm = TRUE)
which(is.na(unlist(significance_year_1_East$significant))) # reference year 1991 is NA, as expected
which(unlist(significance_year_1_East$significant)) + 1990 # 1992-1994 are not significant, 1995-2016 are significant

