# Get CPI and biomass
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(hellometry)


# Load fwdata
fw_database <- 
  readr::read_rds(here::here("data",
                             "scalewebs",
                             "fw_database.rds"))

# Read bromeliad list
bromeliad_list <- 
  readr::read_csv(here::here("data",
                             "bromeliad_list.csv")) %>% 
  dplyr::pull(bromeliad_id)

# Add biomass data --------------------------------------------------------
# Get abundance data from original database
biomass <- 
  fw_database %>% 
  purrr::pluck("abundance") %>% 
  ## Add taxonomy
  hellometry::add_taxonomy() %>% 
  ## Remove columns that are not that important yet
  dplyr::select(-dataset_id, -species_id, -domain, -kingdom,
                -phylum, -subphylum) %>% 
  ## We want dry estimates
  dplyr::mutate(biomass_type = "dry") %>% 
  ## Make column with unknown size class, size and that to receive biomass estimates
  dplyr::mutate(biomass_col = NA,
                size_col = "unknown",
                stage = ifelse(class == "Insecta",
                               "larva", "adult")) %>% 
  ## Add data from Dézerald et al. (2017, full ref above) to have estimates from Cyphon
  ## I used the mean weight from 10-20 individuals with the mean instar size 
  dplyr::bind_rows(tibble::tibble(bwg_name = rep("cypholivier", 10),
                                  class = rep("Insecta", 10),
                                  subclass = rep("Neoptera", 10),
                                  ord = rep("Coleoptera", 10),
                                  subord = rep("Polyphaga", 10),
                                  family = rep("Scirtidae", 10),
                                  genus = rep("Cyphon", 10),
                                  stage = rep("larva", 10),
                                  biomass_type = rep("dry", 10),
                                  biomass_col = c(0.00195, 0.0047, 0.0105, 0.0259, 
                                                  0.0598, 0.0975, 0.1644, 0.3386,
                                                  0.4649, 0.7296), ### NB denominator number of larvae in batch mass measurement
                                  size_col = c("0.14", "0.25", "0.34", "0.44", "0.54",
                                               "0.64", "0.74", "0.84", "0.94", "1.08"),
                                  abundance = c(451, 458, 557, 549, 448, 
                                                430, 338, 263, 221, 58))) %>% 
  ## Do the same for ostracods
  dplyr::bind_rows(tibble::tibble(bwg_name = rep("ostrolivier", 10),
                                  class = rep("Ostracoda", 10),
                                  subclass = rep("Podocopa", 10),
                                  ord = rep("Podocopida", 10),
                                  subord = rep("Cytherocopina", 10),
                                  family = rep("Limnocytheridae", 10),
                                  subfamily = rep("Timiriaseviinae", 10),
                                  genus = rep("Elpidium", 10),
                                  species = rep("Bromeliarum", 10),
                                  stage = rep("adult", 10),
                                  biomass_type = rep("dry", 10),
                                  biomass_col = c(0.0000543, 0.000096914, 0.000180958, 0.000335231, 
                                                  0.000714286, 0.0012, 0.002, 0.003363636,
                                                  0.0055, 0.0172), ### NB denominator number of larvae in batch mass measurement
                                  size_col = c("0.15", "0.24", "0.34", "0.44", "0.54",
                                               "0.64", "0.74", "0.84", "0.94", "1.02"),
                                  abundance = c(18, 141, 447, 765, 1010, 
                                                1779, 3189, 2108, 210, 8)))

# Suite of levels to be used for estimations
level_vec  <-
  c("bwg_name", "genus", "tribe", "subfamily", "family", "subord","ord", 
    "subclass", "class")

# Get the estimations
biomass <-
  hellometry(
    dats = biomass , 
    level_vec  = level_vec , 
    biomass_type = "dry",  
    use_BWG_db = TRUE, 
    no_BWG_data = FALSE,
    r_square_cutoff_upper = 1) ## model for ostracods has a lot of clear data

# Prepare biomass data --------------------------------------------------------
# Basic clean-up 
biomass_clean <- 
  biomass$data %>% 
  dplyr::filter(abundance > 0) %>% 
  dplyr::select(-biomass_type) %>% 
  ## Mass should be in gram
  dplyr::mutate(biomass_col = as.numeric(biomass_col),
                biomass_g = biomass_col/1000) %>% 
  ## Some characters in biomass estimation, remove
  dplyr::filter(!is.na(biomass_col)) %>% 
  ## Remove all the taxonomy and estimation columns
  dplyr::select(-class, -subclass, -ord, -subord, -family,
                -subfamily, -tribe, -genus, -species, -biomass_lower, -biomass_upper,
                -stage, -size_category, -size_level, -size_taxon_name, -model,
                -model_level, -model_taxon_name, -biomass_col) %>% 
  ## Add visit ID
  dplyr::left_join(fw_database %>% 
                     purrr::pluck("bromeliads") %>% 
                     dplyr::select(visit_id, bromeliad_id) %>% 
                     dplyr::mutate(bromeliad_id = as.numeric(bromeliad_id),
                                   visit_id = as.numeric(visit_id)),
                   by = "bromeliad_id") %>% 
  ## Add site codes
  dplyr::left_join(codes %>% 
                     dplyr::select(-replicate, -location, -bio01) %>% 
                     dplyr::distinct(),
                   by = "visit_id") %>%
  ## Remove location
  dplyr::select(-visit_id) %>% 
  ## Sum Coleoptera.48 and Coleoptera.15
  dplyr::group_by(bwg_name, bromeliad_id, site_code) %>% 
  dplyr::summarise(dplyr::across(c(abundance, biomass_g),
                                 ~sum(.))) %>% 
  dplyr::ungroup() %>% 
  ## Get average body mass
  dplyr::mutate(avg_bodymass_g = biomass_g/abundance) %>% 
  ## Rename biomass with total
  dplyr::rename(tot_biomass_g = biomass_g)

# Make sure each row is unique combination of bromeliad id and bwg name
stopifnot(biomass_clean %>% 
            dplyr::select(bromeliad_id, bwg_name) %>% 
            unique() %>% 
            nrow() == biomass_clean %>%
            dplyr::select(bromeliad_id, bwg_name) %>% 
            nrow())

# Add basal resource data, always 1e-11 (lower than lowest biomass),so no real preference
basal_biomass <- 
  biomass_clean %>% 
  ## Select columns and add site code
  dplyr::select(site_code, bromeliad_id) %>% 
  dplyr::distinct() %>% 
  ## Add resource columns
  dplyr::mutate(algae = 1e-11,
                microorganisms = 1e-11,
                FPOM = 1e-11,
                CPOM = 1e-11,
                microzooplankton = 1e-11,
                mesozooplankton = 1e-11,
                bromeliad.tissue = 1e-11) %>% 
  ## Make longer
  tidyr::pivot_longer(cols = c(algae, microorganisms, FPOM,
                               CPOM, microzooplankton, mesozooplankton, bromeliad.tissue),
                      names_to = "bwg_name",
                      values_to = "tot_biomass_g") %>% 
  ## Add 1 as abundance and get avg_bodymass_g
  dplyr::mutate(avg_bodymass_g = tot_biomass_g,
                abundance = 1)

## Combine the two biomass data frames together
biomass_clean <- 
  biomass_clean %>% 
  dplyr::bind_rows(basal_biomass) %>% 
  dplyr::filter(bromeliad_id %in% bromeliad_list)

## Save data
readr::write_csv(biomass_clean,
                 here::here("data",
                            "biomass_clean.csv"))