## CC script

# Libraries
library(tidyverse)
library(here)
library(brms)
library(ggeffects)

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

# Hypothesis 1 - more energy flow with low CPI/high P ------------------------------------------------------------
# Linear model of RELATIVE energy flow entering consumer node ~ CPI and P score
## Model
fluxprop_model_cpip <-
  brms::brm(log(flux_prop + 0.001) ~
              log(cons_CPI_crisp) + log(cons_P) +(1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Save
readr::write_rds(x = fluxprop_model_cpip ,
                 file = here::here("data",
                                   "fluxprop_model_cpip.RDS"))

# Linear model of raw energy flow entering consumer node ~ CPI score
## Model
fluxraw_model_cpip <- 
  brms::brm(log(flux_raw + 0.001) ~
              log(cons_CPI_crisp) + log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
# Save
readr::write_rds(x = fluxraw_model_cpip,
                 file = here::here("data",
                                   "fluxraw_model_cpip.RDS"))

# Hypothesis 2 - more stoichiometric dissimilarity in a pair means less energy flow ------------------------------------------------------------
# With CN distance - relative
fluxprop_model_cndist <-
  brms::brm(log(flux_prop_pair + 0.001) ~
              log(CN_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Save
readr::write_rds(x = fluxprop_model_cndist,
                 file = here::here("data",
                                   "fluxprop_model_cndist.RDS"))

# With CP distance - relative
fluxprop_model_cpdist <- 
  brms::brm(log(flux_prop_pair + 0.001) ~
              log(CP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Save
readr::write_rds(x = fluxprop_model_cpdist,
                 file = here::here("data",
                                   "fluxprop_model_cpdist.RDS"))

# With CNP distance - relative
fluxprop_model_cnpdist <- 
  brms::brm(log(flux_prop_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Save
readr::write_rds(x = fluxprop_model_cnpdist,
                 file = here::here("data",
                                   "fluxprop_model_cnpdist.RDS"))

# With CN distance - raw
fluxraw_model_cndist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CN_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Save
readr::write_rds(x = fluxraw_model_cndist,
                 file = here::here("data",
                                   "fluxraw_model_cndist.RDS"))

# With CP distance - raw
fluxraw_model_cpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Save
readr::write_rds(x = fluxraw_model_cpdist,
                 file = here::here("data",
                                   "fluxraw_model_cpdist.RDS"))

# With CNP distance - raw
fluxraw_model_cnpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
# Save
readr::write_rds(x = fluxraw_model_cnpdist,
                 file = here::here("data",
                                   "fluxraw_model_cnpdist.RDS"))

# Hypothesis 3 - flow if less stoichiometric distances in whole food web --------------------------------
# Total flow CN distance
totalflow_model_cndist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Save
readr::write_rds(x = totalflow_model_cndist,
                 file = here::here("data",
                                   "totalflow_model_cndist.RDS"))

# Total flow CP distance
totalflow_model_cpdist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Save
readr::write_rds(x = totalflow_model_cpdist,
                 file = here::here("data",
                                   "totalflow_model_cpdist.RDS"))

# Total flow CNP distance
totalflow_model_cnpdist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Save
readr::write_rds(x = totalflow_model_cnpdist,
                 file = here::here("data",
                                   "totalflow_model_cnpdist.RDS"))


# Predator flow CN distance
predflow_model_cndist <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Save
readr::write_rds(x = predflow_model_cndist,
                 file = here::here("data",
                                   "predflow_model_cndist.RDS"))

# Predator flow CP distance
predflow_model_cpdist <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Save
readr::write_rds(x = predflow_model_cpdist,
                 file = here::here("data",
                                   "predflow_model_cpdist.RDS"))
# Predator flow CNP distance
predflow_model_cnpdist <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
# Save
readr::write_rds(x = predflow_model_cnpdist,
                 file = here::here("data",
                                   "predflow_model_cnpdist.RDS"))
