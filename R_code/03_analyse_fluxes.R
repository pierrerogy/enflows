# Analysing the fluxes

# Load libraries 
library(igraph)
library(fluxweb)
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(bayestestR)
library(DescTools)
library(ggeffects)

# Read RDS if needed
cpi_model_P <- 
  readr::read_rds(here::here("data",
                             "cpi_model_P.RDS"))
fluxprop_model_p <- 
  readr::read_rds(here::here("data",
                             "fluxprop_model_p.RDS"))
fluxprop_model_cpi <- 
  readr::read_rds(here::here("data",
                             "fluxprop_model_cpi.RDS"))
fluxraw_model_cpi <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_cpi.RDS"))
fluxraw_model_p <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_p.RDS"))
fluxprop_model_cndist <- 
  readr::read_rds(here::here("data",
                             "fluxprop_model_cndist.RDS"))
fluxprop_model_cnpdist <- 
  readr::read_rds(here::here("data",
                             "fluxprop_model_cnpdist.RDS"))
fluxraw_model_cndist <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_cndist.RDS"))
fluxraw_model_cnpdist <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_cnpdist.RDS"))

# Hypothesis 1 - more energy flow with low CPI/high P (P only where measured)------------------------------------------------------------
# Relation between P and CPI
## Histograms
hist(log(fluxes_long_percons$cons_CPI_crisp))
hist(fluxes_long_percons$cons_CPI_crisp)
hist(log(fluxes_long_percons$cons_P))
hist(log(fluxes_long_percons$flux_prop))
hist(log(fluxes_long_percons$flux_raw + 0.01))

## Model 
cpi_model_P <- 
  brms::brm(log(cons_CPI_crisp) ~
              log(avg_bodymass_g) + log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = fluxes_long_percons)
## Look at the effect
bayestestR::describe_posterior(cpi_model_P)

# Linear model of RELATIVE energy flow entering consumer node ~ CPI score
## Model
fluxprop_model_cpi <- 
  brms::brm(log(flux_prop + 0.001) ~
              log(cons_CPI_crisp) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_cpi)

# Linear model of raw energy flow entering consumer node ~ CPI score
## Model
fluxraw_model_cpi <- 
  brms::brm(log(flux_raw + 0.001) ~
              log(cons_CPI_crisp) + (1|bromeliad_id/site),
            iter = 50000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Look at the effect
bayestestR::describe_posterior(fluxraw_model_cpi)

# Linear model of RELATIVE energy flow entering consumer node ~ P content
## Model
fluxprop_model_p <- 
  brms::brm(log(flux_prop + 0.001) ~
              log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = fluxes_long_percons)
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_p)

# Linear model of raw energy flow entering consumer node ~ P content
## Model
fluxraw_model_p <-
  brms::brm(log(flux_raw + 0.001) ~
              log(cons_P) + (1|bromeliad_id/site),
            iter = 50000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Look at the effect
bayestestR::describe_posterior(fluxraw_model_p)
 
# Hypothesis 2 - more stoichiometric dissimilarity in a pair means less energy flow ------------------------------------------------------------
# Histogram
hist(log(fluxes_long_perpair$flux_prop_pair + 0.001))
hist(log(fluxes_long_perpair$flux_raw_pair + 0.001))

# With CN distance - relative
## Model
fluxprop_model_cndist <- 
  brms::brm(log(flux_prop_pair + 0.001) ~
              log(CN_dist) + (1|bromeliad_id/site),
            iter = 50000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_cndist)

# With CNP distance - relative
## Model
fluxprop_model_cnpdist <- 
  brms::brm(log(flux_prop_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 50000,
            family = gaussian(),   
            chains = 3,
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_cnpdist)

# With CN distance - raw
## Model
fluxraw_model_cndist <- 
  brms::brm(log(flux_raw_pair + 0.01) ~
              log(CN_dist) + (1|bromeliad_id/site),
            iter = 50000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Look at the effect
bayestestR::describe_posterior(fluxraw_model_cndist)


# With CNP distance - raw
## Model
fluxraw_model_cnpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 50000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Look at the effect
bayestestR::describe_posterior(fluxraw_model_cnpdist)

# Hypothesis 3 - flow if less stoichiometric distances in whole food web --------------------------------
## Histograms
hist(log(fluxes_long_perbrom$totalflow_Jy))
hist(log(fluxes_long_perbrom$avg_CN_dist))
hist(log(fluxes_long_perbrom$avg_CNP_dist))

# Total flow CN distance
## Model
totalflow_model_CNdist <- 
  brms::brm(log(totalflow_Jy + 0.001) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_CNdist)

# Total flow CNP distance
## Model
totalflow_model_CNPdist <- 
  brms::brm(log(totalflow_Jy + 0.001) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_CNPdist)

# Predator flow CN distance
## Model
predflow_model_CNdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
bayestestR::describe_posterior(predflow_model_CNdist)

# Predator flow CNP distance
## Model
predflow_model_CNPdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Look at the effect
bayestestR::describe_posterior(predflow_model_CNPdist)






