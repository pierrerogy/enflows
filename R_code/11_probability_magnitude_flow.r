# Comparing probability of interaction with energy flow

# Load libraries
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))

# Load BWG database
bwg <- 
        readr::read_rds(here("data", 
                             "scalewebs",
                             "fw_database.rds"))

# Load magnitude of fluxes
fluxes_long <-
    readr::read_csv(here("data",
                         "fluxes_long.csv")) %>%
    ## Remove flux proportion columns
    dplyr::select(-contains("prop"))

# Load all interaction probability matrices
prob_long <-
  load_matrices(folder = here::here("data",
                                    "matrices"),
                ending = "_PM.txt") %>%
  ## Convert each matrix to long format
  purrr::map(~ make_mat_long(matrix = .,
                             unit = "prob_interaction")) %>%
  ## Add site code as column
  purrr::imap(~ dplyr::mutate(.x,
                              site_code = .y)) %>%
  ## Flatten list
  purrr::list_rbind()

# Join probability of interaction to fluxes at the site level
fluxes_prob_long <-
  fluxes_long %>%
  dplyr::left_join(prob_long,
                   by = c("site_code", "consumer", "resource"))

# Load PC scores computed in script 4 (one set of axes per species)
pc_scores <-
  readr::read_csv(here::here("data",
                             "fluxes_long_percons.csv")) %>%
  ## Keep species name and its 4 PC scores
  dplyr::select(bwg_name = consumer,
                dplyr::contains("Axis")) %>%
  dplyr::distinct()

# Add the 4 PC scores for each consumer and resource
fluxes_prob_long <-
  fluxes_prob_long %>%
  ## Consumers
  dplyr::left_join(pc_scores %>%
                     dplyr::rename(consumer = bwg_name) %>%
                     ## to not confuse scores for resources and consumers
                     dplyr::rename_with(~stringr::str_c("con_", .),
                                        dplyr::contains("Axis")),
                   by = "consumer") %>%
  ## Resources
  dplyr::left_join(pc_scores %>%
                     dplyr::rename(resource = bwg_name) %>%
                     dplyr::rename_with(~stringr::str_c("res_", .),
                                        dplyr::contains("Axis")),
                   by = "resource")

