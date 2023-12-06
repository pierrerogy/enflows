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

# Load fluxes
fluxes <- 
  load_matrices(folder = here::here("data",
                                    "fluxes"),
                ending = ".txt")

# Load fwdata
fw_database <- 
  readr::read_rds(here::here("data",
                             "scalewebs",
                             "fw_database.rds"))

# Load stoichimetry data
stoichiometry <- 
  readr::read_csv(here::here("data",
                             "stoichiometry.csv"))


# Resilience/stability ----------------------------------------------------
# Estimategrowth rates of basal sppecies
growth.rates = rep(NA, dim(groups.level$mat)[1])
growth.rates[colSums(groups.level$mat) == 0] = 0.5

#
stability_all <- 
  fluxweb::make.stability(web_matrix, 
                          biomasses = biomasses, 
                          losses = met.rates, 
                          efficiencies = efficiencies)

# Add CPI, P content to flux data -----------------------------------------
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
                                                "bromeliad.tissue", "microorganisms")))
        

# Prepare CPI and P data
traits <- 
  fw_database$traits %>% 
  dplyr::select(bwg_name, dplyr::contains(c("family",
                                            "genus",
                                            "CP", ## CPI
                                            "FD", ## food
                                            "FG", ## feeding mode
                                            "functional_group"
                                            ))) %>% 
  ## Make a mock CPI column
  dplyr::rowwise() %>% 
  dplyr::mutate(CPI_crisp = ifelse(is.na(CP1)|is.na(CP2)|is.na(CP3),
                NA, mean(c(rep(10, CP1), rep(40, CP2), rep(70, CP3)))))

# Combine
fluxes_long <- 
  fluxes_long %>% 
  ## Rename everything according to consumer
  purrr::map(~dplyr::left_join(x = .,
                               y = traits %>% 
                                 dplyr::select(dplyr::everything()) %>% 
                                 dplyr::rename_with(~paste0("cons_",.x )) %>% 
                                 dplyr::rename(consumer = cons_bwg_name),
                               by = "consumer")) %>% 
  purrr::list_rbind() %>% 
  ## Split web column in site and bromeliad
  tidyr::separate_wider_delim(cols = web,
                              delim = "_",
                              names = c("site", "bromeliad_id"))

# Models  -----------------------------------------------------------------
# Mock model with CPI score as the respnse variable
hist(log(fluxes_long$cons_CPI_score))
cpi_model <- 
  brms::brm(log(flux_prop+ 0.01) ~
              log(cons_CPI_score) + (1|bromeliad_id/site),
            iter = 2000,
            family = gaussian(),      
            control = list(adapt_delta = 0.90,
                           max_treedepth = 10),
            data = fluxes_long[1:10000,])
brms::posterior_summary(cpi_model)  
bayestestR::describe_posterior(cpi_model)

