# Follow up analyses
# Load libraries 
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(bayestestR)

# Load site matrices
site_matrices <- 
  load_matrices(folder = here::here("data",
                                    "matrices"),
                ending = "_PM.txt") %>% 
  make_binary_matrices(matrix_list = .,
                       cutoff = 0.02)

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


# Compute models in each site ------------------------------------------
# List to store models and model outputs
outputs_fluxraw_cpip <- 
  list()
outputs_fluxraw_cndist <- 
  list()
outputs_fluxraw_cpdist <- 
  list()
outputs_fluxraw_cnpdist <- 
  list()
outputs_totalflow_cndist <- 
  list()
outputs_totalflow_cpdist <- 
  list()
outputs_totalflow_cnpdist <- 
  list()
outputs_predflow_cndist<- 
  list()
outputs_predflow_cpdist<- 
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
  ### Flux raw and CPI and P
  outputs_fluxraw_cpip[[length(outputs_fluxraw_cpip) + 1]]<-
    brms::brm(log(flux_raw + 0.001) ~
                log(cons_CPI_crisp) + log(cons_P) +(1|bromeliad_id),
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

  ### Flux raw and cn
  outputs_fluxraw_cndist[[length(outputs_fluxraw_cndist) + 1]] <-
    brms::brm(log(flux_raw_pair + 0.001) ~
                log(CN_dist)*lower_cons_CN+ (1|bromeliad_id),
              iter = 20000,
              family = gaussian(),
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Flux raw and cp
  outputs_fluxraw_cpdist[[length(outputs_fluxraw_cpdist) + 1]] <-
    brms::brm(log(flux_raw_pair + 0.001) ~
                log(CP_dist)*lower_cons_CP + (1|bromeliad_id),
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
  ### Total flow and CP
  outputs_totalflow_cpdist[[length(outputs_totalflow_cpdist) + 1]] <-
    brms::brm(log(totalflow_Jy) ~
                log(avg_CP_dist),
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
    brms::brm(log(predflow_Jy + 0.001) ~
                log(avg_CN_dist),
              iter = 20000,
              family = gaussian(),
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Flow to predators and CN
  outputs_predflow_cpdist[[length(outputs_predflow_cpdist) + 1]] <-
    brms::brm(log(predflow_Jy+ 0.001) ~
                log(avg_CP_dist),
              iter = 20000,
              family = gaussian(),
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              data = dats)
  ### Flow to predators and CNP
  outputs_predflow_cnpdist[[length(outputs_predflow_cnpdist) + 1]]<-
    brms::brm(log(predflow_Jy+ 0.001) ~
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
  list(outputs_fluxraw_cpip, outputs_fluxraw_cndist, outputs_fluxraw_cpdist, 
       outputs_fluxraw_cnpdist, outputs_totalflow_cndist, outputs_totalflow_cpdist, 
       outputs_totalflow_cnpdist, outputs_predflow_cndist, 
       outputs_predflow_cpdist, outputs_predflow_cnpdist)
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

# Linear model of raw energy flow entering consumer node ~ CPI
## Model
fluxraw_smallmodel_cpip <-
  brms::brm(log(flux_raw + 0.001) ~
              log(cons_CPI_crisp) + log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_percons %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxraw_smallmodel_cpip)

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
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxraw_smallmodel_cndist)

# With CP distance - raw
## Model
fluxraw_smallmodel_cpdist <- 
  brms::brm(log(flux_raw_pair + 0.01) ~
              log(CP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxraw_smallmodel_cpdist)


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
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(totalflow_smallmodel_cndist)

# Total flow CP distance
## Model
totalflow_smallmodel_cpdist <- 
  brms::brm(log(totalflow_Jy + 0.001) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(totalflow_smallmodel_cpdist)

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
bayestestR::describe_posterior(totalflow_smallmodel_cnpdist)

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
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(predflow_smallmodel_cndist)

# Predator flow CN distance
## Model
predflow_smallmodel_cpdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDC", "CARDO", "PTL")))
## Look at the effect
bayestestR::describe_posterior(predflow_smallmodel_cpdist)

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

# Quadratic term hypothesis 2, and comparing to linear version ---------------------------------------------
# With CN distance - raw
## Model
fluxraw_quadraticmodel_cndist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              I(log(CN_dist)^2)*lower_cons_CN +
              log(CN_dist)*lower_cons_CN + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)

# With CP distance - raw
## Model
fluxraw_quadraticmodel_cpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              I(log(CP_dist)^2)*lower_cons_CP + 
              log(CP_dist)*lower_cons_CP + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)

# With CNP distance - raw
## Model
fluxraw_quadraticmodel_cnpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              I(log(CNP_dist)^2) + log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)

# Compare with LOO
brms::loo(fluxraw_model_cndist, fluxraw_quadraticmodel_cndist)
brms::loo(fluxraw_model_cpdist, fluxraw_quadraticmodel_cpdist)
brms::loo(fluxraw_model_cnpdist, fluxraw_quadraticmodel_cnpdist)

# Quadratic term hypothesis 3 ---------------------------------------------
# Total flow CN distance
## Model
totalflow_quadraticmodel_cndist <- 
  brms::brm(log(totalflow_Jy) ~
              I(log(avg_CN_dist)^2) + log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Total flow CP distance
## Model
totalflow_quadraticmodel_cpdist <- 
  brms::brm(log(totalflow_Jy) ~
              I(log(avg_CP_dist)^2) + log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Total flow CNP distance
## Model
totalflow_quadraticmodel_cnpdist <- 
  brms::brm(log(totalflow_Jy) ~
              I(log(avg_CNP_dist)^2) + log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Predator flow CN distance
## Model
predflow_quadraticmodel_cndist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              I(log(avg_CN_dist)^2) + log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Predator flow CP distance
## Model
predflow_quadraticmodel_cpdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              I(log(avg_CP_dist)^2) + log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Predator flow CNP distance
## Model
predflow_quadraticmodel_cnpdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              I(log(avg_CNP_dist)) + log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)

# Compare with LOO
brms::loo(totalflow_model_cndist, totalflow_quadraticmodel_cndist)
brms::loo(totalflow_model_cpdist, totalflow_quadraticmodel_cpdist)
brms::loo(totalflow_model_cnpdist, totalflow_quadraticmodel_cnpdist)
brms::loo(predflow_model_cndist, predflow_quadraticmodel_cndist)
brms::loo(predflow_model_cpdist, predflow_quadraticmodel_cpdist)
brms::loo(predflow_model_cnpdist, predflow_quadraticmodel_cnpdist)



# Low distance outliers removed, hypothesis 2-------------------
# With CN distance - raw
## Model
fluxraw_nooutliersmodel_cndist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CN_dist)*lower_cons_CN + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(CN_dist >0.001))
## Check assumptions
plot(fluxraw_nooutliersmodel_cndist)
## Look at the effect
bayestestR::describe_posterior(fluxraw_nooutliersmodel_cndist)


# With CP distance - raw
## Model
fluxraw_model_cpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CP_dist)*lower_cons_CP + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(CP_dist >0.001))
## Check assumptions
plot(fluxraw_nooutliersmodel_cpdist)
## Look at the effect
bayestestR::describe_posterior(fluxraw_nooutliersmodel_cpdist)

# With CNP distance - raw
## Model
fluxraw_model_cnpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>% 
              dplyr::filter(CNP_dist >0.001))
## Check assumptions
plot(fluxraw_nooutliersmodel_cnpdist)
## Look at the effect
bayestestR::describe_posterior(fluxraw_nooutliersmodel_cnpdist)


# Low distance outliers removed, hypothesis 3 -----------------------------------
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
## Check assumptions
plot(totalflow_model_cndist_nozero)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cndist_nozero)

# Total flow CN distance
## Model
totalflow_model_cpdist_nozero <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CP_dist > 0.001))
## Check assumptions
plot(totalflow_model_cpdist_nozero)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cpdist_nozero)

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
## Check assumptions
plot(totalflow_model_cnpdist_nozero)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cnpdist_nozero)

# Predator flow CN distance
## Model
predflow_model_cndist_nozero <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CN_dist > 0.001))
## Check assumptions
plot(predflow_model_cndist_nozero)
## Look at the effect
bayestestR::describe_posterior(predflow_model_cndist_nozero)

# Predator flow CP distance
## Model
predflow_model_cpdist_nozero <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CP_dist > 0.001))
## Check assumptions
plot(predflow_model_cpdist_nozero)
## Look at the effect
bayestestR::describe_posterior(predflow_model_cpdist_nozero)

# Predator flow CNP distance
## Model
predflow_model_cnpdist_nozero <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(avg_CNP_dist > 0.001))
## Check assumptions
plot(predflow_model_cnpdist_nozero)
## Look at the effect
bayestestR::describe_posterior(predflow_model_cnpdist_nozero)









