# Analysing the fluxes

# Load libraries 
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(bayestestR)
library(performance)

# Run the 4 MCMC chains of each model in parallel, one per allocated core
# (reads the cores given by SLURM, defaults to 4 for interactive runs)
options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                         unset = "4")))

# Load data
fluxes_long_percons <- 
  readr::read_csv(here::here("data",
                             "fluxes_long_percons.csv"))
fluxes_long_perpair <-
  readr::read_csv(here::here("data",
                             "fluxes_long_perpair.csv"))  
fluxes_long_perbrom <-
  readr::read_csv(here::here("data",
                             "fluxes_long_perbrom.csv"))

# Load models from the cloud ----------------------------------------------
# List of models
 model_list <- 
   c(
     "consflow_model_axis1",
     "consflow_model_axis2",
     "consflow_model_axis3",
     "consflow_model_axis4",
     "pairflow_model_overall",
     "pairflow_model_axis1",
     "pairflow_model_axis2",
     "pairflow_model_axis3",
     "pairflow_model_axis4",
     "totalflow_model_vol",
     "gini_model_vol",
     "totalflow_model_dist",
     "gini_model_dist"
     )

# #  Check if model present in environment, and if not load it
for (model in model_list) {
  if (!exists(model)) {
    assign(model, 
           readr::read_rds(here::here("data",
                                      paste0(model, ".RDS"))))}}


# Hypothesis 1 - models at the scale of consumers ------------------------------------------------------------
## Histograms
hist(fluxes_long_percons$Axis.1)
hist(fluxes_long_percons$Axis.2)
hist(fluxes_long_percons$Axis.3)
hist(fluxes_long_percons$Axis.4)
hist(log(fluxes_long_percons$flux_raw + 1e-3))


# Non-linear model of raw energy flow entering consumer node ~ Axis.1
## Model
consflow_model_axis1 <- 
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.1, 
                    scale = FALSE, center = TRUE) + 
              I(scale(Axis.1, scale = FALSE, center = TRUE)^2) +
              (1 + scale(Axis.1,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Check assumptions
plot(consflow_model_axis1)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis1)
performance::r2_bayes(consflow_model_axis1)

# Non-linear model of raw energy flow entering consumer node ~ Axis.2
## Model
consflow_model_axis2 <- 
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.2, 
                    scale = FALSE, center = TRUE) + 
             I(scale(Axis.2, scale = FALSE, center = TRUE)^2) +
             (1 + scale(Axis.2,
                       scale = FALSE, center = TRUE)|site_code) +
            (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Check assumptions
plot(consflow_model_axis2)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis2)
performance::r2_bayes(consflow_model_axis2)

# Non-linear model of raw energy flow entering consumer node ~ Axis.3
## Model
consflow_model_axis3 <- 
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.3, 
                    scale = FALSE, center = TRUE) + 
             I(scale(Axis.3, scale = FALSE, center = TRUE)^2) +
             (1 + scale(Axis.3,
                       scale = FALSE, center = TRUE)|site_code) +
            (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Check assumptions
plot(consflow_model_axis3)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis3)
performance::r2_bayes(consflow_model_axis3)

# Non-linear model of raw energy flow entering consumer node ~ Axis.4
## Model
consflow_model_axis4 <- 
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.4, 
                    scale = FALSE, center = TRUE) + 
             I(scale(Axis.4, scale = FALSE, center = TRUE)^2) +
             (1 + scale(Axis.4,
                       scale = FALSE, center = TRUE)|site_code) +
            (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Check assumptions
plot(consflow_model_axis4)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis4)
performance::r2_bayes(consflow_model_axis4)


# Save models
readr::write_rds(x = consflow_model_axis1,
                 file = here::here("data",
                                   "consflow_model_axis1.RDS"))
readr::write_rds(x = consflow_model_axis2,
                 file = here::here("data",
                                   "consflow_model_axis2.RDS"))
readr::write_rds(x = consflow_model_axis3,
                 file = here::here("data",
                                   "consflow_model_axis3.RDS"))
readr::write_rds(x = consflow_model_axis4,
                 file = here::here("data",
                                   "consflow_model_axis4.RDS"))

# Hypothesis 2 - more functional dissimilarity in a pair means less energy flow ------------------------------------------------------------
# Histogram
hist(log(fluxes_long_perpair$flux_raw_pair + 0.001))
hist(fluxes_long_perpair$func_dist)
hist(fluxes_long_perpair$dist_axis1)
hist(fluxes_long_perpair$dist_axis2)
hist(fluxes_long_perpair$dist_axis3)
hist(fluxes_long_perpair$dist_axis4)

# Overall functional distance
pairflow_model_overall <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(func_dist,
                    center = TRUE, scale = FALSE) + 
              I(scale(func_dist, center = TRUE, scale = FALSE)^2) +
              (1 + scale(func_dist,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(pairflow_model_overall)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_overall)
performance::r2_bayes(pairflow_model_overall)

# Distance for Axis 1
pairflow_model_axis1 <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis1,
                    center = TRUE, scale = FALSE) + 
            I(scale(dist_axis1, center = TRUE, scale = FALSE)^2) +
            (1 + scale(dist_axis1,
                       center = TRUE, scale = FALSE)|site_code) +
            (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(pairflow_model_axis1)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis1)
performance::r2_bayes(pairflow_model_axis1)

# Distance for Axis 2
pairflow_model_axis2 <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis2,
                    center = TRUE, scale = FALSE) + 
            I(scale(dist_axis2, center = TRUE, scale = FALSE)^2) +
            (1 + scale(dist_axis2,
                       center = TRUE, scale = FALSE)|site_code) +
            (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(pairflow_model_axis2)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis2)
performance::r2_bayes(pairflow_model_axis2)
# Distance for Axis 3
pairflow_model_axis3 <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis3,
                    center = TRUE, scale = FALSE) + 
            I(scale(dist_axis3, center = TRUE, scale = FALSE)^2) +
            (1 + scale(dist_axis3,
                       center = TRUE, scale = FALSE)|site_code) +
            (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(pairflow_model_axis3)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis3)
performance::r2_bayes(pairflow_model_axis3)
# Distance for Axis 4
pairflow_model_axis4 <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis4,
                    center = TRUE, scale = FALSE) + 
            I(scale(dist_axis4, center = TRUE, scale = FALSE)^2) +
            (1 + scale(dist_axis4,
                       center = TRUE, scale = FALSE)|site_code) +
            (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(pairflow_model_axis4)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis4)
performance::r2_bayes(pairflow_model_axis4)

# Save model
readr::write_rds(x = pairflow_model_overall,
                 file = here::here("data",
                                   "pairflow_model_overall.RDS"))
readr::write_rds(x = pairflow_model_axis1,
                 file = here::here("data",
                                   "pairflow_model_axis1.RDS"))
readr::write_rds(x = pairflow_model_axis2,
                 file = here::here("data",
                                   "pairflow_model_axis2.RDS"))
readr::write_rds(x = pairflow_model_axis3,
                 file = here::here("data",
                                   "pairflow_model_axis3.RDS"))
readr::write_rds(x = pairflow_model_axis4,
                 file = here::here("data",
                                   "pairflow_model_axis4.RDS"))


# Hypothesis 3 - more and less imbalanced flow if less functional distances in whole food web --------------------------------
# Histograms
hist(log(fluxes_long_perbrom$totalflow_Jy))
hist(fluxes_long_perbrom$gini)
hist(log(fluxes_long_perbrom$prop_hull))
hist(fluxes_long_perbrom$avg_func_dist)

# Total flow
## Distance
totalflow_model_dist <- 
  brms::brm(log(totalflow_Jy) ~
              scale(avg_func_dist,
                    center = TRUE, scale = FALSE) +
            I(scale(avg_func_dist, center = TRUE, scale = FALSE)^2) +
            (1 + scale(avg_func_dist,
                      center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(totalflow_model_dist)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_dist)
performance::r2_bayes(totalflow_model_dist)
## Volume
totalflow_model_vol <- 
  brms::brm(log(totalflow_Jy) ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
              I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(totalflow_model_vol)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_vol)
performance::r2_bayes(totalflow_model_vol)

# Inequality in flow
## Distance
gini_model_dist <- 
  brms::brm(gini ~
              scale(log(avg_func_dist),
                    center = TRUE, scale = FALSE) +
              I(scale(log(avg_func_dist), center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(avg_func_dist),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(gini_model_dist)
## Look at the effect
bayestestR::describe_posterior(gini_model_dist)
performance::r2_bayes(gini_model_dist)
## Volume
gini_model_vol <- 
  brms::brm(gini ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
               I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2) +
               (1 + scale(log(prop_hull),
                          center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(gini_model_vol)
## Look at the effect
bayestestR::describe_posterior(gini_model_vol)
performance::r2_bayes(gini_model_vol)

# Save models
readr::write_rds(x = totalflow_model_dist,
                 file = here::here("data",
                                   "totalflow_model_dist.RDS"))
readr::write_rds(x = totalflow_model_vol,
                 file = here::here("data",
                                   "totalflow_model_vol.RDS"))
readr::write_rds(x = gini_model_dist,
                 file = here::here("data",
                                   "gini_model_dist.RDS"))
readr::write_rds(x = gini_model_vol,
                 file = here::here("data",
                                   "gini_model_vol.RDS"))




