# Libraries
library(tidyverse)
library(here)
library(brms)
library(ggeffects)

# Load data
# fluxes_long_percons <- 
#   readr::read_csv(here::here("data",
#                              "fluxes_long_percons.csv"))
# fluxes_long_perpair <-
#   readr::read_csv(here::here("data",
#                              "fluxes_long_perpair.csv"))
fluxes_long_perbrom <-
  readr::read_csv(here::here("data",
                             "fluxes_long_perbrom.csv"))


# Models just in sites with stoichiometry ---------------------------------

# # Linear model of RELATIVE energy flow entering consumer node ~ CPI 
# ## Model
# fluxprop_smallmodel_cpi <- 
#   brms::brm(log(flux_prop + 0.001) ~
#               log(cons_CPI_crisp) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 10),
#             data = fluxes_long_percons %>% 
#               dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
# ## Look at the effect
# #bayestestR::describe_posterior(fluxprop_smallmodel_cpi)
# # Save
# readr::write_rds(x = fluxprop_smallmodel_cpi,
#                  file = here::here("data",
#                                    "fluxprop_smallmodel_cpi.RDS"))
# 
# # Linear model of raw energy flow entering consumer node ~ CPI
# ## Model
# fluxraw_smallmodel_cpi <-
#   brms::brm(log(flux_raw + 0.001) ~
#               log(cons_CPI_crisp) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_percons %>% 
#               dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
# ## Look at the effect
# #bayestestR::describe_posterior(fluxraw_smallmodel_cpi)
# # Save
# readr::write_rds(x = fluxprop_smallmodel_cpi,
#                  file = here::here("data",
#                                    "fluxraw_smallmodel_cpi.RDS"))
# 
# # With CN distance - relative
# ## Model
# fluxprop_smallmodel_cndist <- 
#   brms::brm(log(flux_prop_pair + 0.001) ~
#               log(CN_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair %>% 
#               dplyr::filter(site %in% c("CARDC", "CARDO", "PTL", "MAC",
#                                         "SONHI", "SONLOW", "EVDW", "EVPC", 
#                                         "EVTA", "SAUT", "SINN")))
# ## Look at the effect
# #bayestestR::describe_posterior(fluxprop_smallmodel_cndist)
# # Save
# readr::write_rds(x = fluxprop_smallmodel_cpi ,
#                  file = here::here("data",
#                                    "fluxprop_smallmodel_cndist .RDS"))
# 
# # With CNP distance - relative
# ## Model
# fluxprop_smallmodel_cnpdist <- 
#   brms::brm(log(flux_prop_pair + 0.001) ~
#               log(CNP_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),   
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair %>% 
#               dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
# ## Look at the effect
# #bayestestR::describe_posterior(fluxprop_smallmodel_cnpdist)
# # Save
# readr::write_rds(x = fluxprop_smallmodel_cpi ,
#                  file = here::here("data",
#                                    "fluxprop_smallmodel_cnpdist.RDS"))
# 
# # With CN distance - raw
# ## Model
# fluxraw_smallmodel_cndist <- 
#   brms::brm(log(flux_raw_pair + 0.001) ~
#               log(CN_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair %>% 
#               dplyr::filter(site %in% c("CARDC", "CARDO", "PTL", "MAC",
#                                         "SONHI", "SONLOW", "EVDW", "EVPC", 
#                                         "EVTA", "SAUT", "SINN")))
# ## Look at the effect
# #bayestestR::describe_posterior(fluxraw_smallmodel_cndist)
# # Save
# readr::write_rds(x = fluxprop_smallmodel_cpi ,
#                  file = here::here("data",
#                                    "fluxraw_smallmodel_cndist.RDS"))
# 
# 
# # With CNP distance - raw
# ## Model
# fluxraw_smallmodel_cnpdist <- 
#   brms::brm(log(flux_raw_pair + 0.001) ~
#               log(CNP_dist) + (1|bromeliad_id/site),
#             iter = 20000,
#             family = gaussian(),      
#             control = list(adapt_delta = 0.95,
#                            max_treedepth = 15),
#             data = fluxes_long_perpair %>% 
#               dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
# ## Look at the effect
# #bayestestR::describe_posterior(fluxraw_smallmodel_cnpdist)
# # Save
# readr::write_rds(x = fluxprop_smallmodel_cpi ,
#                  file = here::here("data",
#                                    "fluxraw_smallmodel_cnpdist.RDS"))

# Total flow CN distance
## Model
totalflow_smallmodel_cndist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL", "MAC",
                                        "SONHI", "SONLOW", "EVDW", "EVPC", 
                                        "EVTA", "SAUT", "SINN")))
## Look at the effect
#bayestestR::describe_posterior(totalflow_smallmodel_cndist)
# Save
readr::write_rds(x = totalflow_smallmodel_cndist,
                 file = here::here("data",
                                   "totalflow_smallmodel_cndist.RDS"))

# Total flow CNP distance
## Model
totalflow_smallmodel_cnpdist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
#bayestestR::describe_posterior(totalflow_smallmodel_cnpist)
# Save
readr::write_rds(x = totalflow_smallmodel_cnpdist,
                 file = here::here("data",
                                   "totalflow_smallmodel_cnpdist.RDS"))

# Predator flow CN distance
## Model
predflow_smallmodel_cndist <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL", "MAC",
                                        "SONHI", "SONLOW", "EVDW", "EVPC", 
                                        "EVTA", "SAUT", "SINN")))
## Look at the effect
#bayestestR::describe_posterior(predflow_smallmodel_cnist)
# Save
readr::write_rds(x = predflow_smallmodel_cndist,
                 file = here::here("data",
                                   "predflow_smallmodel_cndist.RDS"))

# Predator flow CNP distance
## Model
predflow_smallmodel_cnpdist <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
#bayestestR::describe_posterior(predflow_smallmodel_cnpdist)
# Save
readr::write_rds(x = predflow_smallmodel_cnpdist,
                 file = here::here("data",
                                   "predflow_smallmodel_cnpdist.RDS"))

