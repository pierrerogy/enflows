## CC script

# Libraries
library(tidyverse)
library(here)
library(brms)
library(ggeffects)

# Load data
# fluxes_long_percons <- 
#   readr::read_csv(here::here("data",
#                              "fluxes_long_percons_01.csv"))
# fluxes_long_perpair <- 
#   readr::read_csv(here::here("data",
#                              "fluxes_long_perpair_01.csv"))
fluxes_long_perbrom <- 
  readr::read_csv(here::here("data",
                             "fluxes_long_perbrom_01.csv"))

# # Hypothesis 1 - more energy flow with low CPI/high P (P only where measured)------------------------------------------------------------
# # Linear model of RELATIVE energy flow entering consumer node ~ CPI score
# ## Model
# fluxprop_model_cpi <- 
#   brms::brm(log(flux_prop + 0.001) ~
#               log(cons_CPI_crisp) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.99,
#                            max_treedepth = 15),
#             data = fluxes_long_percons)
# ## Look at the effect
# # bayestestR::describe_posterior(fluxprop_model_cpi)
# # Save
# readr::write_rds(x = fluxprop_model_cpi ,
#                  file = here::here("data",
#                                    "fluxprop_model_cpi_01.RDS"))
# 
# 
# # Linear model of raw energy flow entering consumer node ~ CPI score
# ## Model
# fluxraw_model_cpi <- 
#   brms::brm(log(flux_raw + 0.001) ~
#               log(cons_CPI_crisp) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),   
#             control = list(adapt_delta = 0.99,
#                            max_treedepth = 15),
#             data = fluxes_long_percons)
# ## Look at the effect
# # bayestestR::describe_posterior(fluxraw_model_cpi)
# # Save
# readr::write_rds(x = fluxraw_model_cpi,
#                  file = here::here("data",
#                                    "fluxraw_model_cpi_01.RDS"))
# 
# 
# # Linear model of RELATIVE energy flow entering consumer node ~ P content
# ## Model
# fluxprop_model_p <- 
#   brms::brm(log(flux_prop + 0.001) ~
#               log(cons_P) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 10),
#             data = fluxes_long_percons)
# ## Look at the effect
# # bayestestR::describe_posterior(fluxprop_model_p)
# # Save
# readr::write_rds(x = fluxprop_model_p,
#                  file = here::here("data",
#                                    "fluxprop_model_p_01.RDS"))
# 
# # Linear model of raw energy flow entering consumer node ~ P content
# ## Model
# fluxraw_model_p <-
#   brms::brm(log(flux_raw + 0.001) ~
#               log(cons_P) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_percons)
# ## Look at the effect
# bayestestR::describe_posterior(fluxraw_model_p)
# # Save
# readr::write_rds(x = fluxraw_model_p,
#                  file = here::here("data",
#                                    "fluxraw_model_p_01.RDS"))
# # Hypothesis 2 - more stoichiometric dissimilarity in a pair means less energy flow ------------------------------------------------------------
# # With CN distance - relative
# ## Model
# fluxprop_model_cndist <- 
#   brms::brm(log(flux_prop_pair + 0.001) ~
#               log(CN_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair)
# ## Look at the effect
# # bayestestR::describe_posterior(fluxprop_model_cndist)
# # Save
# readr::write_rds(x = fluxprop_model_cndist,
#                  file = here::here("data",
#                                    "fluxprop_model_cndist_01.RDS"))
# 
# # With CNP distance - relative
# ## Model
# fluxprop_model_cnpdist <- 
#   brms::brm(log(flux_prop_pair + 0.001) ~
#               log(CNP_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),   
#             chains = 3,
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair)
# ## Look at the effect
# # bayestestR::describe_posterior(fluxprop_model_cnpdist)
# # Save
# readr::write_rds(x = fluxprop_model_cnpdist,
#                  file = here::here("data",
#                                    "fluxprop_model_cnpdist_01.RDS"))
# 
# # With CN distance - raw
# ## Model
# fluxraw_model_cndist <- 
#   brms::brm(log(flux_raw_pair + 0.01) ~
#               log(CN_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair)
# ## Look at the effect
# # bayestestR::describe_posterior(fluxraw_model_cndist)
# # Save
# readr::write_rds(x = fluxraw_model_cndist,
#                  file = here::here("data",
#                                    "fluxraw_model_cndist_01.RDS"))
# 
# 
# # With CNP distance - raw
# ## Model
# fluxraw_model_cnpdist <- 
#   brms::brm(log(flux_raw_pair + 0.001) ~
#               log(CNP_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair)
# ## Look at the effect
# # bayestestR::describe_posterior(fluxraw_model_cnpdist)
# # Save
# readr::write_rds(x = fluxraw_model_cnpdist,
#                  file = here::here("data",
#                                    "fluxraw_model_cnpdist_01.RDS"))
# 
# Hypothesis 3 - flow if less stoichiometric distances in whole food web --------------------------------
# Total flow CN distance
## Model
totalflow_model_cndist <- 
  brms::brm(log(totalflow_Jy + 0.001) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
# bayestestR::describe_posterior(totalflow_model_CNdist)
# Save
readr::write_rds(x = totalflow_model_cndist,
                 file = here::here("data",
                                   "totalflow_model_cndist_01.RDS"))

# Total flow CNP distance
## Model
totalflow_model_cnpdist <- 
  brms::brm(log(totalflow_Jy + 0.001) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
# bayestestR::describe_posterior(totalflow_model_cnpdist)
# Save
readr::write_rds(x = totalflow_model_cnpdist,
                 file = here::here("data",
                                   "totalflow_model_cnpdist_01.RDS"))


# Predator flow CN distance
## Model
predflow_model_cndist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
# bayestestR::describe_posterior(predflow_model_CNdist)
# Save
readr::write_rds(x = predflow_model_cndist,
                 file = here::here("data",
                                   "predflow_model_cndist_01.RDS"))

# Predator flow CNP distance
## Model
predflow_model_cnpdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
# bayestestR::describe_posterior(predflow_model_CNPdist)
# Save
readr::write_rds(x = predflow_model_cnpdist,
                 file = here::here("data",
                                   "predflow_model_cnpdist_01.RDS"))
