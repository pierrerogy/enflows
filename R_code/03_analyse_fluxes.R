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

# Read RDS if needed
## Hypothesis 1
cpi_model_P <- 
  readr::read_rds(here::here("data",
                             "cpi_model_P.RDS"))
fluxraw_model_cpip <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_cpip.RDS"))
## Hypothesis 2
fluxraw_model_cndist <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_cndist.RDS"))
fluxraw_model_cpdist <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_cpdist.RDS"))
fluxraw_model_cnpdist <- 
  readr::read_rds(here::here("data",
                             "fluxraw_model_cnpdist.RDS"))
## Hypothesis 3
totalflow_model_cndist <- 
  readr::read_rds(here::here("data",
                             "totalflow_model_cndist.RDS"))
totalflow_model_cpdist <- 
  readr::read_rds(here::here("data",
                             "totalflow_model_cpdist.RDS"))
totalflow_model_cnpdist <- 
  readr::read_rds(here::here("data",
                             "totalflow_model_cnpdist.RDS"))
predflow_model_cndist <- 
  readr::read_rds(here::here("data", 
                             "predflow_model_cndist.RDS"))
predflow_model_cpdist <- 
  readr::read_rds(here::here("data", 
                             "predflow_model_cpdist.RDS"))
predflow_model_cnpdist <- 
  readr::read_rds(here::here("data",
                             "predflow_model_cnpdist.RDS"))

# Hypothesis 1 - more energy flow with low CPI/high P ------------------------------------------------------------
# Relation between P and CPI
## Histograms
hist(log(fluxes_long_percons$cons_CPI_crisp))
hist(log(fluxes_long_percons$cons_P))
hist(log(fluxes_long_percons$flux_raw + 1e-3))

## Model 
cpi_model_P <- 
  brms::brm(log(cons_CPI_crisp) ~
              log(avg_bodymass_g) + log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = fluxes_long_percons)
## Check assumptions
plot(cpi_model_P)
## Look at the effect
bayestestR::describe_posterior(cpi_model_P)
brms::bayes_R2(cpi_model_P)

# Linear model of raw energy flow entering consumer node ~ CPI score and P
## Model
fluxraw_model_cpip <- 
  brms::brm(log(flux_raw + 0.001) ~
              log(cons_CPI_crisp) + log(cons_P) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),   
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Check assumptions
plot(fluxraw_model_cpip)
## Look at the effect
bayestestR::describe_posterior(fluxraw_model_cpip)
brms::bayes_R2(fluxraw_model_cpip)
 
# Hypothesis 2 - more stoichiometric dissimilarity in a pair means less energy flow ------------------------------------------------------------
# Histogram
hist(log(fluxes_long_perpair$flux_raw_pair + 0.001))

# With CN distance - raw
## Model
fluxraw_model_cndist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CN_dist)*lower_cons_CN + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(fluxraw_model_cndist)
## Look at the effect
### Get contrasts
con <- 
  ## Create reference grid to be worked on by emmeans
  emmeans::ref_grid(fluxraw_model_cndist) %>% 
  ## Get pairwise differences
  emmeans:::pairs.emmGrid() 
bayestestR::describe_posterior(con)
brms::bayes_R2(fluxraw_model_cndist)


# With CP distance - raw
## Model
fluxraw_model_cpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CP_dist)*lower_cons_CP + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(fluxraw_model_cpdist)
## Look at the effect
bayestestR::describe_posterior(fluxraw_model_cpdist)
brms::bayes_R2(fluxraw_model_cpdist)

# With CNP distance - raw
## Model
fluxraw_model_cnpdist <- 
  brms::brm(log(flux_raw_pair + 0.001) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Check assumptions
plot(fluxraw_model_cnpdist)
## Look at the effect
bayestestR::describe_posterior(fluxraw_model_cnpdist)
brms::bayes_R2(fluxraw_model_cnpdist)


# Hypothesis 3 - flow if less stoichiometric distances in whole food web --------------------------------
## Histograms
hist(log(fluxes_long_perbrom$totalflow_Jy))
hist(log(fluxes_long_perbrom$predflow_Jy))
hist(log(fluxes_long_perbrom$avg_CN_dist))
hist(log(fluxes_long_perbrom$avg_CNP_dist))

# Total flow CN distance
## Model
totalflow_model_cndist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(totalflow_model_cndist)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cndist)
brms::bayes_R2(totalflow_model_cndist)


# Total flow CP distance
## Model
totalflow_model_cpdist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(totalflow_model_cpdist)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cpdist)
brms::bayes_R2(totalflow_model_cpdist)

# Total flow CNP distance
## Model
totalflow_model_cnpdist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(totalflow_model_cnpdist)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_cnpdist)
brms::bayes_R2(totalflow_model_cnpdist)

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
## Check assumptions
plot(predflow_model_cndist)
## Look at the effect
bayestestR::describe_posterior(predflow_model_cndist)
brms::bayes_R2(predflow_model_cndist)

# Predator flow CP distance
## Model
predflow_model_cpdist <- 
  brms::brm(log(predflow_Jy + 0.001) ~
              log(avg_CP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom)
## Check assumptions
plot(predflow_model_cpdist)
## Look at the effect
bayestestR::describe_posterior(predflow_model_cpdist)
brms::bayes_R2(predflow_model_cpdist)

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
## Check assumptions
plot(predflow_model_cnpdist)
## Look at the effect
brms::bayes_R2(predflow_model_cnpdist)



