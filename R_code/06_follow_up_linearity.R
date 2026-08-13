# Follow up analyses - testing whether the relationships in script 05 are
# better described by a linear term.

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

# Load non-linear models from script 05 ---------------------------------------
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

# Hypothesis 1 - linear models at the scale of consumers --------------
# Same as script 05 but with a linear fixed effect.
# Linear flow ~ Axis.1
consflow_model_axis1_linear <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.1, scale = FALSE, center = TRUE) +
              (1 + scale(Axis.1,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Linear flow ~ Axis.2
consflow_model_axis2_linear <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.2, scale = FALSE, center = TRUE) +
              (1 + scale(Axis.2,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Linear flow ~ Axis.3
consflow_model_axis3_linear <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.3, scale = FALSE, center = TRUE) +
              (1 + scale(Axis.3,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Linear flow ~ Axis.4
consflow_model_axis4_linear <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.4, scale = FALSE, center = TRUE) +
              (1 + scale(Axis.4,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)

# Compare each non-linear model with its linear equivalent (same response)
brms::loo(consflow_model_axis1, 
          consflow_model_axis1_linear)
brms::loo(consflow_model_axis2, 
          consflow_model_axis2_linear)
brms::loo(consflow_model_axis3, 
          consflow_model_axis3_linear)
brms::loo(consflow_model_axis4, 
          consflow_model_axis4_linear)

# Save models and remove from memorey to make room
## Save
readr::write_rds(consflow_model_axis1_linear,
                 here::here("data",
                            "consflow_model_axis1_linear.RDS"))
readr::write_rds(consflow_model_axis2_linear,
                 here::here("data",
                            "consflow_model_axis2_linear.RDS"))
readr::write_rds(consflow_model_axis3_linear,
                 here::here("data",
                            "consflow_model_axis3_linear.RDS"))
readr::write_rds(consflow_model_axis4_linear,
                 here::here("data",
                            "consflow_model_axis4_linear.RDS"))
## Delete
rm(consflow_model_axis1_linear,
   consflow_model_axis1,
   consflow_model_axis2_linear,
   consflow_model_axis2,
   consflow_model_axis3_linear,
   consflow_model_axis3, 
   consflow_model_axis4_linear,
   consflow_model_axis4)

# Hypothesis 2 - linear models at the scale of consumer-resource pairs ----
# Linear flow ~ overall functional distance
pairflow_model_overall_linear <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(func_dist, center = TRUE, scale = FALSE) +
              (1 + scale(func_dist,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Linear flow ~ distance for Axis 1
pairflow_model_axis1_linear <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis1, center = TRUE, scale = FALSE) +
              (1 + scale(dist_axis1,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Linear flow ~ distance for Axis 2
pairflow_model_axis2_linear <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis2, center = TRUE, scale = FALSE) +
              (1 + scale(dist_axis2,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Linear flow ~ distance for Axis 3
pairflow_model_axis3_linear <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis3, center = TRUE, scale = FALSE) +
              (1 + scale(dist_axis3,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Linear flow ~ distance for Axis 4
pairflow_model_axis4_linear <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis4, center = TRUE, scale = FALSE) +
              (1 + scale(dist_axis4,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)

# Compare each non-linear model with its linear equivalent (same response)
brms::loo(pairflow_model_overall, 
          pairflow_model_overall_linear)
brms::loo(pairflow_model_axis1, 
          pairflow_model_axis1_linear)
brms::loo(pairflow_model_axis2, 
          pairflow_model_axis2_linear)
brms::loo(pairflow_model_axis3, 
          pairflow_model_axis3_linear)
brms::loo(pairflow_model_axis4, 
          pairflow_model_axis4_linear)

# Save models and remove from memorey to make room
## Save
readr::write_rds(pairflow_model_overall_linear,
                 here::here("data",
                            "pairflow_model_overall_linear.RDS"))
readr::write_rds(pairflow_model_axis1_linear,
                 here::here("data",
                            "pairflow_model_axis1_linear.RDS"))
readr::write_rds(pairflow_model_axis2_linear,
                 here::here("data",
                            "pairflow_model_axis2_linear.RDS"))
readr::write_rds(pairflow_model_axis3_linear,
                 here::here("data",
                            "pairflow_model_axis3_linear.RDS"))
readr::write_rds(pairflow_model_axis4_linear,
                 here::here("data",
                            "pairflow_model_axis4_linear.RDS"))
## Delete
rm(pairflow_model_overall_linear, 
   pairflow_model_overall,
   pairflow_model_axis1_linear, 
   pairflow_model_axis1,
   pairflow_model_axis2_linear, 
   pairflow_model_axis2,
   pairflow_model_axis3_linear,
   pairflow_model_axis3,
   pairflow_model_axis4_linear
   pairflow_model_axis4)

# Hypothesis 3 - linear models at the scale of whole food webs --------
# Total flow - volume
totalflow_model_vol_linear <-
  brms::brm(log(totalflow_Jy) ~
              scale(log(prop_hull), center = TRUE, scale = FALSE) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Total flow - distance
totalflow_model_dist_linear <-
  brms::brm(log(totalflow_Jy) ~
              scale(avg_func_dist, center = TRUE, scale = FALSE) +
              (1 + scale(avg_func_dist,
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Flow to predators - volume
predflow_model_vol_linear <-
  brms::brm(log(predflow_Jy + 0.01) ~
              scale(log(prop_hull), center = TRUE, scale = FALSE) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Flow to predators - distance
predflow_model_dist_linear <-
  brms::brm(log(predflow_Jy + 0.01) ~
              scale(avg_func_dist, center = TRUE, scale = FALSE) +
              (1 + scale(avg_func_dist,
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Gini index - volume
gini_model_vol_linear <-
  brms::brm(gini ~
              scale(log(prop_hull), center = TRUE, scale = FALSE) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Gini index - distance
gini_model_dist_linear <-
  brms::brm(gini ~
              scale(log(avg_func_dist), center = TRUE, scale = FALSE) +
              (1 + scale(log(avg_func_dist),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Compare each non-linear model with its linear equivalent (same response)
brms::loo(totalflow_model_vol, 
          totalflow_model_vol_linear)
brms::loo(totalflow_model_dist, 
          totalflow_model_dist_linear)
brms::loo(predflow_model_vol, 
          predflow_model_vol_linear)
brms::loo(predflow_model_dist, 
          predflow_model_dist_linear)
brms::loo(gini_model_vol, 
          gini_model_vol_linear)
brms::loo(gini_model_dist, 
          gini_model_dist_linear)


# Save the linear models
readr::write_rds(totalflow_model_vol_linear,
                 here::here("data",
                            "totalflow_model_vol_linear.RDS"))
readr::write_rds(totalflow_model_dist_linear,
                 here::here("data",
                            "totalflow_model_dist_linear.RDS"))
readr::write_rds(predflow_model_vol_linear,
                 here::here("data",
                            "predflow_model_vol_linear.RDS"))
readr::write_rds(predflow_model_dist_linear,
                 here::here("data",
                            "predflow_model_dist_linear.RDS"))
readr::write_rds(gini_model_vol_linear,
                 here::here("data",
                            "gini_model_vol_linear.RDS"))
readr::write_rds(gini_model_dist_linear,
                 here::here("data",
                            "gini_model_dist_linear.RDS"))