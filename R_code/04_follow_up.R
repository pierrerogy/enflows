# Follow up analyses
# Load libraries 
library(igraph)
library(fluxweb)
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(bayestestR)
library(cowplot)


# Compute models in each site ------------------------------------------
# List to store models and model outputs
outputs_fluxprop_CPI <- 
  list()
outputs_fluxprop_p <- 
  list()
outputs_fluxraw_CPI <- 
  list()
outputs_fluxraw_p <- 
  list()
outputs_fluxprop_cndist <- 
  list()
outputs_fluxprop_cnpdist <- 
  list()
outputs_fluxraw_cndist <- 
  list()
outputs_fluxraw_cnpdist <- 
  list()
outputs_totalflow_cndist <- 
  list()
outputs_totalflow_cnpdist <- 
  list()
outputs_predflow_cndist<- 
  list()
outputs_predflow_cnpdist <- 
  list()

# Loop through each site to get models
for(i in names(site_matrices)){
  ## Sanity check
  print(i)
  ## First do models per cons
  dats <- 
    fluxes_long_percons %>% 
    dplyr::filter(site == i)
  ### Flux prop and CPI
  outputs_fluxprop_CPI[[length(outputs_fluxprop_CPI) + 1]] <- 
    brms::brm(log(flux_prop + 0.001) ~
                log(cons_CPI_crisp) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Flux raw and CPI
  outputs_fluxraw_CPI[[length(outputs_fluxraw_CPI) + 1]]<- 
    brms::brm(log(flux_raw + 0.001) ~
                log(cons_CPI_crisp) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),   
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Sanity check
  print(i)
  ### Flux prop and P
  outputs_fluxprop_p[[length(outputs_fluxprop_p) + 1]]  <- 
    brms::brm(log(flux_prop + 0.001) ~
                log(cons_P) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 10),
              data = dats)
  ### Flux raw and P
  outputs_fluxraw_p[[length(outputs_fluxraw_p) + 1]]  <-    
    brms::brm(log(flux_raw + 0.001) ~
                log(cons_P) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  
  ## Sanity check 
  print(i)
  
  ## Second do models per pair
  dats <- 
    fluxes_long_perpair %>% 
    dplyr::filter(site == i)
  ### Flux prop and cn distance
  outputs_fluxprop_cndist[[length(outputs_fluxprop_cndist) + 1]] <- 
    brms::brm(log(flux_prop_pair + 0.001) ~
                log(CN_dist) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Flux prop and cnp distance
  outputs_fluxprop_cnpdist[[length(outputs_fluxprop_cnpdist) + 1]]  <- 
    brms::brm(log(flux_prop_pair + 0.001) ~
                log(CNP_dist) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),   
              chains = 3,
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Sanity check
  print(i)
  ### Flux raw and cn 
  outputs_fluxraw_cndist[[length(outputs_fluxraw_cndist) + 1]] <- 
    brms::brm(log(flux_raw_pair + 0.001) ~
                log(CN_dist) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Flux raw and CNP
  outputs_fluxraw_cnpdist[[length(outputs_fluxraw_cnpdist) + 1]] <- 
    brms::brm(log(flux_raw_pair + 0.001) ~
                log(CNP_dist) + (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  
  ## Sanity check
  print(i)
  
  ## Third do models per bromeliad
  dats <- 
    fluxes_long_perbrom %>% 
    dplyr::filter(site == i)
  ### Total flow and CN
  outputs_totalflow_cndist[[length(outputs_totalflow_cndist) + 1]] <- 
    brms::brm(log(totalflow_Jy) ~
                log(avg_CN_dist),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Total flow and CNP
  outputs_totalflow_cnpdist[[length(outputs_totalflow_cnpdist) + 1]]  <- 
    brms::brm(log(totalflow_Jy) ~
                log(avg_CNP_dist),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  
  ### Sanity check
  print(i)
  ### Flow to predators and CN
  outputs_predflow_cndist[[length(outputs_predflow_cndist) + 1]] <- 
    brms::brm(log(predflow_Jy) ~
                log(avg_CN_dist),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Flow to predators and CNP
  outputs_predflow_cnpdist[[length(outputs_predflow_cnpdist) + 1]]<- 
    brms::brm(log(predflow_Jy) ~
                log(avg_CNP_dist),
              iter = 20000,
              family = gaussian(),      
              control = list(adapt_delta = 0.95,
                             max_treedepth = 15),
              data = dats)
}

# Save all lists as single list
## Make
masterlist_sites <- 
  list(outputs_fluxprop_CPI, outputs_fluxprop_p, outputs_fluxraw_CPI, outputs_fluxraw_p,
       outputs_fluxprop_cndist, outputs_fluxprop_cnpdist, outputs_fluxraw_cndist, outputs_fluxraw_cnpdist,
       outputs_totalflow_cndist, outputs_totalflow_cnpdist, outputs_predflow_cndist, outputs_predflow_cnpdist)
## Save
readr::write_rds(masterlist_sites,
                 here::here("data",
                            "masterlist_sites.RDS"))


# Models just in sites with stoichiometry ---------------------------------
## Model 
cpi_smallmodel_P <- 
  brms::brm(log(cons_CPI_crisp) ~
              log(avg_bodymass_g) + log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = fluxes_long_percons %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(cpi_smallmodel_P)

# Linear model of RELATIVE energy flow entering consumer node ~ P content
## Model
fluxprop_smallmodel_p <- 
  brms::brm(log(flux_prop + 0.001) ~
              log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = fluxes_long_percons %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxprop_smallmodel_p)

# Linear model of raw energy flow entering consumer node ~ P content
## Model
fluxraw_smallmodel_p <-
  brms::brm(log(flux_raw + 0.001) ~
              log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_percons %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxraw_smallmodel_p)

# Linear model of RELATIVE energy flow entering consumer node ~ CPI 
## Model
fluxprop_smallmodel_cpi <- 
  brms::brm(log(flux_prop + 0.001) ~
              log(cons_CPI_crisp) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 10),
            data = fluxes_long_percons %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxprop_smallmodel_cpi)

# Linear model of raw energy flow entering consumer node ~ CPI
## Model
fluxraw_smallmodel_cpi <-
  brms::brm(log(flux_raw + 0.001) ~
              log(cons_CPI_crisp) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_percons %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxraw_smallmodel_cpi)

# With CN distance - relative
## Model
fluxprop_smallmodel_cndist <- 
  brms::brm(log(flux_prop_pair + 0.001) ~
              log(CN_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL", "MAC",
                                        "SONHI", "SONLOW", "EVDW", "EVPC", 
                                        "EVTA", "SAUT", "SINN")))
## Look at the effect
bayestestR::describe_posterior(fluxprop_smallmodel_cndist)

# With CNP distance - relative
## Model
fluxprop_smallmodel_cnpdist <- 
  brms::brm(log(flux_prop_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxprop_smallmodel_cnpdist)

# With CN distance - raw
## Model
fluxraw_smallmodel_cndist <- 
  brms::brm(log(flux_raw_pair + 0.01) ~
              log(CN_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL", "MAC",
                                        "SONHI", "SONLOW", "EVDW", "EVPC", 
                                        "EVTA", "SAUT", "SINN")))
## Look at the effect
bayestestR::describe_posterior(fluxraw_smallmodel_cndist)


# With CNP distance - raw
## Model
fluxraw_smallmodel_cnpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxraw_smallmodel_cnpdist)

# Total flow CN distance
## Model
totalflow_smallmodel_cndist <- 
  brms::brm(log(totalflow_Jy + 0.001) ~
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
bayestestR::describe_posterior(totalflow_smallmodel_cndist)

# Total flow CNP distance
## Model
totalflow_smallmodel_cnpdist <- 
  brms::brm(log(totalflow_Jy + 0.001) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(totalflow_smallmodel_cnpist)

# Predator flow CN distance
## Model
predflow_smallmodel_cndist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
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
bayestestR::describe_posterior(predflow_smallmodel_cnist)

# Predator flow CNP distance
## Model
predflow_smallmodel_cnpdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(predflow_smallmodel_cnpdist)

# Bromeliad-scale models but lowest values removed to avoid having that little bump-----------------------------------
# Total flow CN distance
## Model
totalflow_model_cndist_nozero <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CN_dist > 0.001))
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cndist_nozero)

# Total flow CNP distance
## Model
totalflow_model_cnpdist_nozero <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CNP_dist > 0.001)) 
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cnpdist_nozero)

# Predator flow CN distance
## Model
predflow_model_cndist_nozero <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CN_dist > 0.001))

## Look at the effect
bayestestR::describe_posterior(predflow_model_cndist_nozero)

# Predator flow CNP distance
## Model
predflow_model_cnpdist_nozero <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CNP_dist > 0.001))
## Look at the effect
bayestestR::describe_posterior(predflow_model_cnpdist_nozero)



