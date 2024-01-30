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

# Prep data ---------------------------------------------------------------
# Basal species
basal <- 
  c("algae","FPOM","CPOM","microorganisms",
    "microzooplankton","mesozooplankton", "bromeliad.tissue")  

# Load fluxes
fluxes <- 
  load_matrices(folder = here::here("data",
                                    "fluxes"),
                ending = ".txt")

# Load site matrices
site_matrices <- 
  load_matrices(folder = here::here("data",
                                    "matrices"),
                ending = "_PM.txt") %>% 
  make_binary_matrices(matrix_list = .,
                       cutoff = 0.02)

# Load CPI and stoichimetry data
CPI_stoich <-
  readr::read_csv(here::here("data",
                             "CPI_stoich.csv"))

# Convert matrices to list
fluxes_long <- 
  fluxes %>% 
  purrr::map(~ make_mat_long(matrix = .,
                             unit = "flux")) %>% 
  ## Rename flux column to raw
  purrr::map(~ dplyr::rename(.,
                             flux_raw = flux)) %>% 
  ## Convert flux data to proportion for each web
  purrr::map(~ dplyr::mutate(.,
                             flux_prop = flux_raw/sum(flux_raw))) %>% 
  ## Add name of web as column
  purrr::imap(~dplyr::mutate(.x,
                             web = .y)) %>% 
  ## Just remove those rows where consumers are resources
  purrr::map(~ dplyr::filter(.,
                             consumer %notin% c("algae", "FPOM", "CPOM",
                                                "bromeliad.tissue", "microorganisms"))) %>% 
  ## Flatten list
  purrr::list_rbind() %>% 
  ## Split web column in site and bromeliad
  tidyr::separate_wider_delim(cols = web,
                              delim = "_",
                              names = c("site", "bromeliad_id")) %>% 
  ## Remove cannibalism
  dplyr::filter(consumer != resource)
  
# Combine for each consumer
fluxes_long_percons <- 
  fluxes_long %>% 
  ## Pool per consumer
  dplyr::select(-resource) %>% 
  dplyr::group_by(site, bromeliad_id, consumer) %>% 
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~sum(.))) %>%
  dplyr::ungroup() %>% 
  ## Add trait data relative to consumer
  dplyr::left_join(CPI_stoich %>% 
                     dplyr::select(dplyr::everything()) %>% 
                     dplyr::rename_with(~paste0("cons_",.x )) %>% 
                     dplyr::rename(consumer = cons_bwg_name),
                   by = "consumer") %>% 
  ## Add biomass data relative to consumers
  dplyr::left_join(readr::read_csv(here::here("data",
                                              "biomass_clean.csv")) %>% 
                     dplyr::select(consumer = bwg_name, site = site_code ,bromeliad_id,
                                   tot_biomass_g, avg_bodymass_g, abundance) %>% 
                     dplyr::mutate(bromeliad_id = as.character(bromeliad_id)),
                   by = c("site", "bromeliad_id", "consumer")) %>% 
  ## Remove proportion 0 (e.g. small food webs)
  dplyr::filter(flux_prop > 0)

# Get distance between consumer/resource pairs
fluxes_long_perpair <- 
  fluxes_long %>% 
  ## Remove forbidden links (see Srivastava et al. (Scalewebs))
  dplyr::anti_join(
    readr::read_csv(here::here("data",
                               "forbidden_links.csv")) %>% 
      dplyr::filter(is == 0) %>% 
      tidyr::separate(col = interaction_id,
                      into = c("site", "n"),
                      sep = "\\.") %>% 
      dplyr::select(site, consumer, resource),
    by = c("site", "consumer", "resource")) %>% 
  ## Remove interactions where we have basal species
  dplyr::filter(resource %notin% basal & consumer %notin% basal) %>% 
  ## Add the composition of consumer
  dplyr::left_join(CPI_stoich %>% 
                     dplyr::select(bwg_name, functional_group, CPI_crisp,
                                   C_mean, N_mean, P_mean) %>% 
                     dplyr::rename_with(~paste0("cons_",.x )) %>% 
                     dplyr::rename(consumer = cons_bwg_name),
                   by = "consumer") %>% 
  ## Add the composition of resource
  dplyr::left_join(CPI_stoich %>% 
                     dplyr::select(bwg_name, functional_group, CPI_crisp,
                                   C_mean, N_mean, P_mean) %>% 
                     dplyr::rename_with(~paste0("res_",.x )) %>% 
                     dplyr::rename(resource = res_bwg_name),
                   by = "resource") %>% 
  ## Get CN distance
  dplyr::mutate(CN_dist = sqrt((cons_C_mean - res_C_mean)^2 + (cons_N_mean - res_N_mean)^2)) %>% 
  ## Get CNP distance
  dplyr::mutate(CNP_dist = sqrt((cons_C_mean - res_C_mean)^2 + (cons_N_mean - res_N_mean)^2 + (cons_P_mean - res_P_mean)^2)) %>% 
  ## Remove columns
  dplyr::select(-dplyr::contains("mean")) %>% 
  ## Replace distances of 0 with Nas
  dplyr::mutate(CN_dist = ifelse(CN_dist == 0, 
                                 NA, CN_dist),
                CNP_dist = ifelse(CNP_dist == 0, 
                                 NA, CNP_dist))



# Get Gini coefficient and total energy flux for each bromeliad, 
# based on all the fluxes (i.e. not pooled per consumer, excluding zooplankton)
fluxes_long_perbrom <- 
  fluxes_long_percons %>% 
  ## Remove zooplankton
  dplyr::filter(!stringr::str_detect(string = consumer,
                                     pattern = "plankton")) %>% 
  dplyr::select(site, bromeliad_id, flux_raw, cons_C_mean, cons_N_mean, cons_P_mean) %>%
  ## Group by site and bromeliad id
  dplyr::group_by(site, bromeliad_id) %>%
  ## Summarise
  dplyr::summarise(dplyr::across(flux_raw,
                                 ~ DescTools::Gini(.),
                                 .names = "gini"),
                   dplyr::across(flux_raw,
                                 ~ sum(.),
                                 .names = "totalflow_Jy")) %>% 
  dplyr::ungroup() %>% 
  ## Add flow just to predators
  dplyr::left_join(fluxes_long_percons %>% 
                     ## Keep just predators
                     dplyr::filter(cons_functional_group %in% c("piercer", "engulfer")) %>%
                     dplyr::select(site, bromeliad_id, flux_raw) %>%
                     ## Group by site and bromeliad id
                     dplyr::group_by(site, bromeliad_id) %>%
                     ## Summarise
                     dplyr::summarise(dplyr::across(flux_raw,
                                                    ~ sum(.),
                                                    .names = "predflow_Jy")),
                   by = c("site", "bromeliad_id")) %>% 
  ## Add average CN and CNP stoichiometric distance
  dplyr::left_join(fluxes_long_percons %>% 
                     ## Remove zooplankton
                     dplyr::filter(!stringr::str_detect(string = consumer,
                                                        pattern = "plankton")) %>% 
                     dplyr::select(site, bromeliad_id, flux_raw, cons_C_mean, cons_N_mean, cons_P_mean) %>%
                     ## Group by site and bromeliad id
                     dplyr::nest_by(site, bromeliad_id) %>% 
                     ## Get average
                     dplyr::mutate(avg_CN_dist = average_dist(data, CN = TRUE),
                                   avg_CNP_dist = average_dist(data, CNP = TRUE)) %>% 
                     ## Select columns
                     dplyr::select(-data),
                   by = c("site", "bromeliad_id")) %>% 
  ## Remove NaN and 0 
  dplyr::mutate(dplyr::across(c(avg_CN_dist, avg_CNP_dist),
                              ~ifelse(is.nan(.),
                                     NA, ifelse(. == 0,
                                                NA, .)))) 
  

# Hypothesis 1 - more energy flow with low CPI/high P (P only where measured)------------------------------------------------------------
# Relation between P and CPI
## Histograms
hist(log(fluxes_long_percons$cons_CPI_crisp))
hist(fluxes_long_percons$cons_CPI_crisp)
hist(log(fluxes_long_percons %>% 
           dplyr::filter(site %in% c("CARDO", "CARDC", "PTL")) %>% 
           dplyr::select(cons_P_mean) %>%
           dplyr::pull()))
## Model 
cpi_model_P <- 
  brms::brm(log(cons_CPI_crisp) ~
              log(avg_bodymass_g) + log(cons_P_mean) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = fluxes_long_percons %>% 
              dplyr::filter(site %in% c("CARDO", "CARDC", "PTL")))
## Look at the effect
bayestestR::describe_posterior(cpi_model_P)
## Plot by hand because here covariate
# Get effect
model_effect <- 
  brms::conditional_effects(cpi_model_P,
                            method = "fitted")$cons_P_mean %>% 
  dplyr::mutate(estimate__ = exp(estimate__),
                lower__ = exp(lower__),
                upper__ = exp(upper__))

# Linear model of RELATIVE energy flow entering consumer node ~ CPI score
## Histogram
hist(sqrt(fluxes_long_percons$flux_prop))
## Model
fluxprop_model_cpi <- 
  brms::brm(sqrt(flux_prop) ~
              log(cons_CPI_crisp) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_percons)
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_cpi)

# Linear model of RELATIVE energy flow entering consumer node ~ P content
## Model
fluxprop_model_p <- 
  brms::brm(sqrt(flux_prop) ~
              log(cons_P_mean) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.92,
                           max_treedepth = 10),
            data = fluxes_long_percons %>%
              dplyr::filter(site %in% c("CARDO", "CARDC", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_p)
  
# Hypothesis 2 - more stoichiometric dissimilarity in a pair means less energy flow ------------------------------------------------------------
# Histogram
hist(sqrt(fluxes_long_perpair$flux_prop))

# With CN distance
## Model
fluxprop_model_cndist <- 
  brms::brm(sqrt(flux_prop) ~
              log(CN_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair)
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_cndist)


# With CNP distance
## Model
fluxprop_model_cnpdist <- 
  brms::brm(sqrt(flux_prop) ~
              log(CNP_dist) + (1|bromeliad_id/site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair %>%
              dplyr::filter(site %in% c("CARDO", "CARDC", "PTL")))
## Look at the effect
bayestestR::describe_posterior(fluxprop_model_cnpdist)

# Hypothesis 3 - flow if less stoichiometric distances in whole food web --------------------------------
# Detritivore dissimilarity vs total energy flow
## Histograms
hist(log(fluxes_long_perbrom$totalflow_Jy))
hist(log(fluxes_long_perbrom$avg_CN_dist))
hist(log(fluxes_long_perbrom$avg_CNP_dist))

# Total flow CN distance
## Model
totalflow_model_CNdist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(totalflow_Jy > 0))
## Look at the effect
bayestestR::describe_posterior(totalflow_model_CNdist)

# Total flow CNP distance
## Model
totalflow_model_CNPdist <- 
  brms::brm(log(totalflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDO", "CARDC", "PTL") &
                                          totalflow_Jy > 0))
## Look at the effect
bayestestR::describe_posterior(totalflow_model_CNPdist)

# Predator flow CN distance
## Model
predflow_model_CNdist <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CN_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(totalflow_Jy > 0))
## Look at the effect
bayestestR::describe_posterior(predflow_model_CNdist)

# Predator flow CNP distance
## Model
predflow_model_CNPdist <- 
  brms::brm(log(predflow_Jy) ~
              log(avg_CNP_dist) + (1|site),
            iter = 20000,
            family = gaussian(),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom %>% 
              dplyr::filter(site %in% c("CARDO", "CARDC", "PTL") &
                              totalflow_Jy > 0))
## Look at the effect
bayestestR::describe_posterior(predflow_model_CNPdist)






# Hypothesis 3 - polynomial relationship between dissimilarity and productivity  ------------------------------------------------------------
# Hypothesis 3.2 - polynomial relationship between dissimilarity and energy flow to predators  ------------------------------------------------------------
# Detritivore dissimilarity vs total energy flow
## Histograms
hist(log(fluxes_long_perbrom_CABO$predflow_Jy))
hist(log(fluxes_long_perbrom_CABO$detr_dist))

## Model 
predflow_model_detrdist <- 
  brms::brm(log(predflow_Jy) ~
              poly(log(detr_dist), degree = 2) #+ (1|/site)
            ,
            iter = 5000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = fluxes_long_perbrom_CABO)
## Look at the effect
bayestestR::describe_posterior(predflow_model_detrdist)
## Plot
basic_line(model = predflow_model_detrdist, 
           dats = fluxes_long_perbrom_CABO, 
           x_var = "detr_dist", 
           y_var = "totalflow_Jy",
           scale = "log")




## Reflects evenness of energy use by consumer
## community and energetic overdominance. High
## inequality is associated with low biodiversity and
## stability of the system

# Average 3D C,N,P distance between consumers vs. Gini index
## Histograms
hist(log(fluxes_long_perbrom_CABO$gini))
hist(log(fluxes_long_perbrom_CABO$detr_dist))
## Model
gini_model_detrdist <- 
  brms::brm(log(gini) ~
              log(detr_dist) #+ (1|/site)
            ,
            iter = 5000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = fluxes_long_perbrom_CABO)
## Look at the effect
bayestestR::describe_posterior(gini_model_detrdist)
## Plot
basic_line(model = gini_model_detrdist, 
           dats = fluxes_long_perbrom_CABO, 
           x_var = "detr_dist", 
           y_var = "gini",
           scale = "log")


