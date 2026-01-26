# Soilfoodwebs attempt
library(soilfoodwebs)
library(tidyverse)
library(here)

# Look at the data that we need for soilfoodwebs
intro_comm
## Structure
str(intro_comm)
## one predator/prey matrix
## prey are columns, predator rows
intro_comm$imat
## Then one dataframe of 10 columns
intro_comm$prop
### ID: species name
### d: Death rate
### a: Assimilation efficiency, carbon assimilated per carbon consumed
### p: Produciton efficiency, biomass gained per carbon assimilated
### B: Equilibrium biomass
### CN: Fixed carbon to nitrogen ratio of a pool
### DetritusRecycling: The proportion of detritus produced in the food web sent to that node
### isDetritus: detritus or not, binary 0/1
### isPlant: Does the pool gain outside carbon based on a per biomass growth rate?, binary 0/1
### canIMM: Can the node immobilize inorganic nitrogen? Should the function corrstoich ignore it?

# Extract a small food web to play with from script 02
web <- 
  site_matrices[[1]]

speclist <- 
 foodwebs_full %>% 
  dplyr::filter(bromeliad_id == "12181",
                site_code == "CABO") %>% 
  dplyr::select(consumer, resource)
speclist <- 
  unique(c(speclist$consumer, speclist$resource))

## Make into a foodweb
web_1 <- 
  web[speclist, speclist]

## Assimilation efficiency from Benoît's calculations, for CABO site only
temp.kT <- 
  0.5391269
assim <- 
  exp(2.266) * exp(0.164 * temp.kT) / (1 + exp(2.266) * exp(0.164 * temp.kT))

## Get other info
web_2 <- 
  readr::read_csv(here::here("data",
                             "biomass_stoich.csv")) %>% 
  dplyr::select(bwg_name, bromeliad_id,site_code,
                tot_biomass_g,C, N, P) %>% 
  dplyr::filter(bromeliad_id == "12181",
                site_code == "CABO") %>% 
  ## Add functional group
  dplyr::left_join(fw_database$traits %>% 
                      dplyr::select(bwg_name, 
                                    functional_group)) %>%
  ## Add assimilation efficiencies based on trophic level
  dplyr::mutate(a = ifelse(is.na(functional_group),
                           1, 
                           assim)) %>% 
  ## Add nutrient ratios
  dplyr::mutate(CN = C/N,
                CP = C/P,
                NP = N/P) %>%
  ## Add detritus recycling and Boolean
  dplyr::mutate(DetritusRecycling = ifelse(bwg_name == "FPOM",
                                           1,
                                           0),
                isDetritus = ifelse(stringr::str_detect(string = bwg_name,
                                                        pattern = "POM"),
                                    1, 
                                    0),
                isPlant = ifelse(bwg_name == "algae",
                                 1,
                                 0),
                canIMM = ifelse(is.na(functional_group),
                                1, 
                                0)) %>%
  dplyr::select(ID = bwg_name,
                B = tot_biomass_g,
                a,
                CN, CP, NP,
                DetritusRecycling,
                isDetritus,
                isPlant,
                canIMM)
 
  


## DOUBLE CHEACK ASSIMILATION EFFFICIENCY
  
  # death rates
  # prodution efficiency