# Follow up analyses - testing whether the relationships in script 05 are
# better described by a non-linear (second-order polynomial) term.

# Load libraries
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(brms)

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

# Load linear models from script 05 ---------------------------------------
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
    "predflow_model_vol",
    "gini_model_vol",
    "totalflow_model_dist",
    "predflow_model_dist",
    "gini_model_dist"
    )

#  Check if model present in environment, and if not load it
for (model in model_list) {
  if (!exists(model)) {
    assign(model,
           readr::read_rds(here::here("data",
                                       paste0(model, ".RDS"))))}}

# Hypothesis 1 - non-linear models at the scale of consumers --------------
# Same as script 05 but with a quadratic fixed effect (x + I(x^2));
# a raw squared term is used rather than poly() so that brms drops rows with
# missing predictor values (poly() errors on NAs), keeping the LOO comparison
# valid because both models are then fit on the same observations
# Non-linear flow ~ Axis.1
consflow_model_axis1_nl <-
  brms::brm(log(flux_raw + 0.001) ~
              Axis.1 + I(Axis.1^2) +
              (1 + scale(Axis.1,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Non-linear flow ~ Axis.2
consflow_model_axis2_nl <-
  brms::brm(log(flux_raw + 0.001) ~
              Axis.2 + I(Axis.2^2) +
              (1 + scale(Axis.2,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Non-linear flow ~ Axis.3
consflow_model_axis3_nl <-
  brms::brm(log(flux_raw + 0.001) ~
              Axis.3 + I(Axis.3^2) +
              (1 + scale(Axis.3,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Non-linear flow ~ Axis.4
consflow_model_axis4_nl <-
  brms::brm(log(flux_raw + 0.001) ~
              Axis.4 + I(Axis.4^2) +
              (1 + scale(Axis.4,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)

# Compare each linear model with its non-linear equivalent (same response)
brms::loo(consflow_model_axis1, 
          consflow_model_axis1_nl)
brms::loo(consflow_model_axis2, 
          consflow_model_axis2_nl)
brms::loo(consflow_model_axis3, 
          consflow_model_axis3_nl)
brms::loo(consflow_model_axis4, 
          consflow_model_axis4_nl)

# Hypothesis 2 - non-linear models at the scale of consumer-resource pairs ----
# Non-linear flow ~ overall functional distance
pairflow_model_overall_nl <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              func_dist + I(func_dist^2) +
              (1 + scale(func_dist,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Non-linear flow ~ distance for Axis 1
pairflow_model_axis1_nl <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              dist_axis1 + I(dist_axis1^2) +
              (1 + scale(dist_axis1,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Non-linear flow ~ distance for Axis 2
pairflow_model_axis2_nl <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              dist_axis2 + I(dist_axis2^2) +
              (1 + scale(dist_axis2,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Non-linear flow ~ distance for Axis 3
pairflow_model_axis3_nl <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              dist_axis3 + I(dist_axis3^2) +
              (1 + scale(dist_axis3,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Non-linear flow ~ distance for Axis 4
pairflow_model_axis4_nl <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              dist_axis4 + I(dist_axis4^2) +
              (1 + scale(dist_axis4,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)

# Compare each linear model with its non-linear equivalent (same response)
brms::loo(pairflow_model_overall, 
          pairflow_model_overall_nl)
brms::loo(pairflow_model_axis1, 
          pairflow_model_axis1_nl)
brms::loo(pairflow_model_axis2, 
          pairflow_model_axis2_nl)
brms::loo(pairflow_model_axis3, 
          pairflow_model_axis3_nl)
brms::loo(pairflow_model_axis4, 
          pairflow_model_axis4_nl)

# Hypothesis 3 - non-linear models at the scale of whole food webs --------
# Total flow - volume
totalflow_model_vol_nl <-
  brms::brm(log(totalflow_Jy) ~
              log(prop_hull) + I(log(prop_hull)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Total flow - distance
totalflow_model_dist_nl <-
  brms::brm(log(totalflow_Jy) ~
              avg_func_dist + I(avg_func_dist^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Flow to predators - volume
predflow_model_vol_nl <-
  brms::brm(log(predflow_Jy + 0.01) ~
              log(prop_hull) + I(log(prop_hull)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Flow to predators - distance
predflow_model_dist_nl <-
  brms::brm(log(predflow_Jy + 0.01) ~
              avg_func_dist + I(avg_func_dist^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Gini index - volume
gini_model_vol_nl <-
  brms::brm(gini ~
              log(prop_hull) + I(log(prop_hull)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Gini index - distance
gini_model_dist_nl <-
  brms::brm(gini ~
              log(avg_func_dist) + I(log(avg_func_dist)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Compare each linear model with its non-linear equivalent (same response)
brms::loo(totalflow_model_vol, 
          totalflow_model_vol_nl)
brms::loo(totalflow_model_dist, 
          totalflow_model_dist_nl)
brms::loo(predflow_model_vol, 
          predflow_model_vol_nl)
brms::loo(predflow_model_dist, 
          predflow_model_dist_nl)
brms::loo(gini_model_vol, 
          gini_model_vol_nl)
brms::loo(gini_model_dist, 
          gini_model_dist_nl)


# Save the non-linear models
readr::write_rds(consflow_model_axis1_nl,
                 here::here("data",
                            "consflow_model_axis1_nl.RDS"))
readr::write_rds(consflow_model_axis2_nl,
                 here::here("data",
                            "consflow_model_axis2_nl.RDS"))
readr::write_rds(consflow_model_axis3_nl,
                 here::here("data",
                            "consflow_model_axis3_nl.RDS"))
readr::write_rds(consflow_model_axis4_nl,
                 here::here("data",
                            "consflow_model_axis4_nl.RDS"))
readr::write_rds(pairflow_model_overall_nl,
                 here::here("data",
                            "pairflow_model_overall_nl.RDS"))
readr::write_rds(pairflow_model_axis1_nl,
                 here::here("data",
                            "pairflow_model_axis1_nl.RDS"))
readr::write_rds(pairflow_model_axis2_nl,
                 here::here("data",
                            "pairflow_model_axis2_nl.RDS"))
readr::write_rds(pairflow_model_axis3_nl,
                 here::here("data",
                            "pairflow_model_axis3_nl.RDS"))
readr::write_rds(pairflow_model_axis4_nl,
                 here::here("data",
                            "pairflow_model_axis4_nl.RDS"))
readr::write_rds(totalflow_model_vol_nl,
                 here::here("data",
                            "totalflow_model_vol_nl.RDS"))
readr::write_rds(totalflow_model_dist_nl,
                 here::here("data",
                            "totalflow_model_dist_nl.RDS"))
readr::write_rds(predflow_model_vol_nl,
                 here::here("data",
                            "predflow_model_vol_nl.RDS"))
readr::write_rds(predflow_model_dist_nl,
                 here::here("data",
                            "predflow_model_dist_nl.RDS"))
readr::write_rds(gini_model_vol_nl,
                 here::here("data",
                            "gini_model_vol_nl.RDS"))
readr::write_rds(gini_model_dist_nl,
                 here::here("data",
                            "gini_model_dist_nl.RDS"))