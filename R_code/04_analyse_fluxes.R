# Analysing the fluxes

# Load libraries 
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(bayestestR)

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


# Hypothesis 1 - more energy flow to species with rapid life cycles (low CPI) ------------------------------------------------------------
# Relation between P and CPI
## Histograms
hist(log(fluxes_long_percons$cons_CPI_defuzzified))
hist(log(fluxes_long_percons$flux_raw + 1e-3))


# Linear model of raw energy flow entering consumer node ~ CPI score 
## Model
consflow_model <- 
  brms::brm(log(flux_raw + 0.001) ~
              scale(log(cons_CPI_defuzzified), 
                    scale = FALSE, center = TRUE) + 
              (1 + scale(log(cons_CPI_defuzzified),
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Check assumptions
plot(consflow_model)
## Look at the effect
bayestestR::describe_posterior(consflow_model)
brms::bayes_R2(consflow_model)


# Linear model of CPI score ~ average body mass 
## Model
cpi_model <- 
  brms::brm(log(cons_CPI_defuzzified) ~
              scale(log(avg_bodymass_g), 
                    scale = FALSE, center = TRUE),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons %>% 
              dplyr::select(consumer, avg_bodymass_g, cons_CPI_defuzzified) %>% 
              dplyr::distinct())
## Check assumptions
plot(cpi_model)
## Look at the effect
bayestestR::describe_posterior(cpi_model)
brms::bayes_R2(cpi_model)
 
# Save models
readr::write_rds(x = consflow_model,
                 file = here::here("data",
                                   "consflow_model.RDS"))
readr::write_rds(x = cpi_model,
                 file = here::here("data",
                                   "cpi_model.RDS"))
# Hypothesis 2 - more functional dissimilarity in a pair means less energy flow ------------------------------------------------------------
# Histogram
hist(log(fluxes_long_perpair$flux_raw_pair + 0.001))
hist(fluxes_long_perpair$func_dist )


# Model
pairflow_model <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(func_dist,
                    center = TRUE, scale = FALSE)+ 
              (1 + scale(func_dist,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(pairflow_model)
## Look at the effect
bayestestR::describe_posterior(pairflow_model)
brms::bayes_R2(pairflow_model)

# Save model
readr::write_rds(x = pairflow_model,
                 file = here::here("data",
                                   "pairflow_model.RDS"))

# Hypothesis 3 - more and less imbalanced flow if less functional distances in whole food web --------------------------------
# Histograms
hist(log(fluxes_long_perbrom$totalflow_Jy))
hist(log(fluxes_long_perbrom$predflow_Jy))
hist(fluxes_long_perbrom$gini)
hist(log(fluxes_long_perbrom$prop_hull))
hist(fluxes_long_perbrom$avg_func_dist)

# Total flow
## Volume
totalflow_model_vol <- 
  brms::brm(log(totalflow_Jy) ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
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
brms::bayes_R2(totalflow_model_vol)
## Distance
totalflow_model_dist <- 
  brms::brm(log(totalflow_Jy) ~
              scale(avg_func_dist,
                    center = TRUE, scale = FALSE) +
            (1 + scale(log(prop_hull),
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
brms::bayes_R2(totalflow_model_dist)

# Flow to predators
## Volume
predflow_model_vol <- 
  brms::brm(log(predflow_Jy + 0.01) ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
            (1 + scale(log(prop_hull),
                      center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(predflow_model_vol)
## Look at the effect
bayestestR::describe_posterior(predflow_model_vol)
brms::bayes_R2(predflow_model_vol)
## Distance
predflow_model_dist <- 
  brms::brm(log(predflow_Jy + 0.01) ~
              scale(avg_func_dist,
                    center = TRUE, scale = FALSE) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(predflow_model_dist)
## Look at the effect
bayestestR::describe_posterior(predflow_model_dist)
brms::bayes_R2(predflow_model_dist)

# Flow to predators
## Volume
gini_model_vol <- 
  brms::brm(gini ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
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
brms::bayes_R2(gini_model_vol)
## Distance
gini_model_dist <- 
  brms::brm(gini ~
              scale(log(avg_func_dist),
                    center = TRUE, scale = FALSE) +
              (1 + scale(log(prop_hull),
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
brms::bayes_R2(gini_model_dist)

# Save models
readr::write_rds(x = totalflow_model_vol,
                 file = here::here("data",
                                   "totalflow_model_vol.RDS"))
readr::write_rds(x = predflow_model_vol,
                 file = here::here("data",
                                   "predflow_model_vol.RDS"))
readr::write_rds(x = gini_model_vol,
                 file = here::here("data",
                                   "gini_model_vol.RDS"))
readr::write_rds(x = totalflow_model_dist,
                 file = here::here("data",
                                   "totalflow_model_dist.RDS"))
readr::write_rds(x = predflow_model_dist,
                 file = here::here("data",
                                   "predflow_model_dist.RDS"))
readr::write_rds(x = gini_model_dist,
                 file = here::here("data",
                                   "gini_model_dist.RDS"))

