# Calculating energy fluxes in the food webs
## Author: Pierre Rogy (based on a script from Benoît Gauzens)

# Load libraries 
library(igraph)
library(fluxweb)
library(tidyverse)
library(here)
library(hellometry)
source(here::here("R_code",
                  "functions.R"))

# Load site codes
codes <- 
  readr::read_csv(here::here("data", "scalewebs",
                             "site_codes_visits.csv")) %>% 
  dplyr::select(-`...1`)

# Load fwdata
fw_database <- 
  readr::read_rds(here::here("data",
                             "scalewebs",
                             "fw_database.rds"))


# Get biomass data --------------------------------------------------------
# Read data
dats <- 
  readr::read_csv(here::here("data", "scalewebs",
                             "abundance_cleaned.csv")) 

# Add taxonomy to data
dats <- 
  hellometry::add_taxonomy(data_table = dats)

# Add columns required by the package to work
dats <-
  readr::read_csv(here::here("data", "scalewebs",
                             "abundance_cleaned.csv")) %>% 
  hellometry::add_taxonomy() %>% 
  ## Everyone's size will be unknown
  dplyr::mutate(size_mm = "unknown") %>%
  ## Fix the couple problems with missing taxonomy
  ### two coleoptera with that info only in bwg_name
  dplyr::mutate(ord = ifelse(bwg_name %in% c("Coleoptera.19", "Coleoptera.48"),
                             "Coleoptera", ord),
                class = ifelse(bwg_name %in% c("Coleoptera.19", "Coleoptera.48"),
                               ## some lost platyhelminthes
                               "Insecta", ifelse(phylum == "Platyhelminthes",
                                                 "Turbellaria", class))) %>% 
  ## Add column with development stage
  ## Let's make everything all insects larvae and rest adults
  dplyr::mutate(stage = ifelse(class == "Insecta",
                               "larva", "adult"))

# Quick check to make sure things are alright   
str(dats)

# Add biomass
dats_biomass <- 
  hellometry::hello_metry(data_table = dats,
                          biomass_kind = "dry",
                          database = TRUE)

# Basic clean-up 
biomass_clean <- 
  dats_biomass %>% 
  dplyr::filter(abundance > 0) %>% 
  dplyr::select(-size_original, -path, 
                -stage, -size_used_mm,
                -visit_id, -dataset_id) %>% 
  ## Mass should be in gram
  dplyr::mutate(biomass_mg = as.numeric(biomass_mg),
                biomass_g = biomass_mg/1000) %>% 
  ## Some characters in biomass estimation, remove
  dplyr::filter(!is.na(biomass_mg)) %>% 
  ## Remove the biomass_mg and the confidence intervals
  dplyr::select(-biomass_mg, -ci_lwr, -ci_upr) %>% 
  ## Remove all the taxonomy columns
  dplyr::select(-species_id, -domain:-species) %>% 
  ## Sum Coleoptera.48 and Coleoptera.15
  dplyr::group_by(bwg_name, bromeliad_id, location) %>% 
  dplyr::summarise(dplyr::across(c(abundance, biomass_g),
                                 ~sum(.))) %>% 
  dplyr::ungroup() %>% 
  ## Get average body mass
  dplyr::mutate(avg_bodymass_g = biomass_g/abundance) %>% 
  ## Rename biomass with total
  dplyr::rename(tot_biomass_g = biomass_g) %>% 
  ## Add site codes
  dplyr::left_join(codes %>% 
                     dplyr::select(-replicate, -visit_id, -bio01) %>% 
                     dplyr::distinct(),
                   by = "location") %>%
  ## Remove location
  dplyr::select(-location)

# Make sure each row is unique combination of bromeliad id and bwg name
stopifnot(biomass_clean %>% 
            dplyr::select(bromeliad_id, bwg_name) %>% 
            unique() %>% 
            nrow() == biomass_clean %>%
            dplyr::select(bromeliad_id, bwg_name) %>% 
            nrow())

# Add basal resource data, always 3.456e-11 (lower than lowest biomass),so no real preference
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
  dplyr::bind_rows(basal_biomass)

## Save data
readr::write_csv(biomass_clean,
                 here::here("data",
                            "biomass_clean.csv"))

# Get food webs -----------------------------------------------------------
# Read biomass data
biomass_clean <- 
  readr::read_csv(here::here("data",
                             "biomass_clean.csv"))


#  Register all site-level matrices 
site_matrices <- 
  load_matrices(folder = here::here("data",
                                    "matrices"),
                ending = "_PM.txt") %>% 
  make_binary_matrices(matrix_list = .,
                       cutoff = 0.02)

# Common inputs to all matrices
## Basal species
basal <- 
  c("algae","FPOM","CPOM","microorganisms",
    "microzooplankton","mesozooplankton", "bromeliad.tissue")  
## List of taxon names
taxon_bwg <- 
  readr::read_csv(here::here("data", "scalewebs",
                             "taxon_bwg_lookup.csv")) %>% 
  ## add basal species
  dplyr::bind_rows(tibble::tibble(taxon_name = basal,
                                     bwg_name = basal))
  
# Compute food webs
foodwebs <- 
  make_food_webs(taxon_names = taxon_bwg,
                 basal = basal,
                 binary_matrices = site_matrices,
                 biomass_data = biomass_clean)

# Get temperature data ----------------------------------------------------
# Add Boltzmann constant 
boltz <- 
  0.00008617343

# Add temperature at for each species pool, convert to K
temp_data <- 
  codes %>% 
  dplyr::select(site_code, bio01) %>% 
  dplyr::rename(mean_temp = bio01) %>% 
  dplyr::mutate(mean_temp = mean_temp + 273.15) %>% 
  dplyr::distinct()

## Use 20 degrees as standardised temperature across all dataset
T0 <- 
  273.15 + 20   

# Combine data ------------------------------------------------------------

# Combine
foodwebs_full <- 
  foodwebs %>% 
  # Combine food webs and temperature
  dplyr::left_join(temp_data,
                   by = "site_code") %>% 
  # Combine food webs with biomass and body mass data for consumers
  dplyr::left_join(biomass_clean %>% 
                     dplyr::select(-abundance) %>% 
                     dplyr::rename(consumer = bwg_name,
                                   con_bio = tot_biomass_g,
                                   con_mass = avg_bodymass_g),
                   by = c("consumer", "bromeliad_id", "site_code")) %>% 
  # Combine food webs with biomass and body mass data for resources
  dplyr::left_join(biomass_clean %>% 
                     dplyr::select(-abundance) %>% 
                     dplyr::rename(resource = bwg_name,
                                   res_bio = tot_biomass_g,
                                   res_mass = avg_bodymass_g),
                   by = c("resource", "bromeliad_id", "site_code")) %>% 
  ## Remove cannibalism
  dplyr::filter(consumer != resource)



# Compute fluxes - code adapted from Benoît Gauzens ----------------------------------------------------------
# Get list of unique webs (i.e. bromeliads)
webs <- 
  unique(foodwebs_full$bromeliad_id)

# If the folder to save the results does not exist, create it
if (!dir.exists(here::here("data", "fluxes"))){
  dir.create(here::here("data", "fluxes"))
}

# Make a loop over all the different networks
## interaction matrix, vector of efficiency and metabolism, vector of biomass use body mass and T for met
for (web in webs){
  
  # Subset the data based on a specific food web
  tab <- 
    foodwebs_full %>% 
    dplyr::filter(bromeliad_id == web)
  # Get temperature at site
  temp  <- 
    tab %>% 
    dplyr::select(mean_temp) %>% 
    dplyr::distinct() %>% 
    dplyr::pull()
  # Convert temperature to kT  
  temp.kT <- 
    (temp - T0) /(boltz * (temp) * T0)
  # Get site name
  site <- 
    tab %>% 
    dplyr::select(site_code) %>% 
    dplyr::distinct() %>% 
    dplyr::pull()
  
  # As the network is based on columns of preys and predators I use the functions from igraph to 
  # create the network matrix
  ## Get list of species
  species <- 
    tab %>% 
    dplyr::select(species = resource) %>% 
    dplyr::bind_rows(tab %>% 
                       dplyr::select(species = consumer)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull()
    
  ## Make list of edges and convert to network and matrix
  web_matrix <- 
    tab %>% 
    dplyr::select(resource, consumer) %>% 
    as.matrix() %>% 
    ### Conversion using igraph
    igraph::graph_from_edgelist(directed = T) %>% 
    igraph::as_adjacency_matrix(sparse = FALSE)
  ### Reorder the matrix accordingly to the order in species
  ### !!!!!! important
  order <- 
    match(species, 
          colnames(web_matrix))
  web_matrix <- 
    web_matrix[order, order] 
  
  # Now get the body masses and biomasses of species (essentially redoing the biomass_clean data)
  list.bms <- 
    tab %>% 
    dplyr::select(species = consumer,
                  BM = con_mass,
                  biomass = con_bio) %>% 
    dplyr::bind_rows(tab %>% 
                       dplyr::select(species = resource,
                                     BM = res_mass,
                                     biomass = res_bio)) %>% 
    dplyr::distinct()
    
  # Check if each species has unique values
  ## Gives error if not
  stopifnot(length(list.bms$species) == length(unique(list.bms$species)))

  # Reorder the array accordingly to the vector species
  list.bms <- 
    list.bms[match(species, list.bms$species), ]
  
  # Deal with basal resources
  ## Detritus and bromeliad tissue
  dets <- 
    c("FPOM", "CPOM")
  pl <- 
    c("bromeliad.tissue", "algae")
  ## Save the detritus and plant tissue that are actually present in the web
  detritus <- 
    species[species %in% dets]
  plants <- 
    species[species %in% pl]
  ## Put the two inside sources
  sources <- 
    species[species %in% plants | species %in% detritus]
  ## Microorganims are included in animals
  animals <- 
    species[(!(species %in% dets)) & (!(species %in% plants))]
  ## Extract non basal species, ie those that do not consume anything 
  non.basals <-  
    rowSums(web_matrix) > 0.0
  
  # Extract all biomasses
  biomasses  <- 
    list.bms %>% 
    dplyr::select(biomass) %>% 
    ## set the biomass of species for which there is no information to 1
    ## this is not important for the amount of flux in the web, but to 
    ## define the preferences of species on these resources
    dplyr::mutate(biomass = ifelse(is.na(biomass), 
                                           1, biomass)) %>% 
    ## Convert to vector
    dplyr::pull()
  
  # Get metabolic rate of species without body mass to 1, this is not important 
  # if and only if these species are basal species
  met.rates  <- 
    exp((-0.29 * log(list.bms$BM) + 17.17) - 0.69 / (boltz * (273.15 + temp)))
  met.rates[is.na(met.rates)] <-  
    1 
  
  # Set all assimilatin efficiencies as for animal prey
  efficiencies <-  
    rep(exp(2.266) * exp(0.164 * temp.kT) / (1 + exp(2.266) * exp(0.164 * temp.kT)), 
        length(species))
  # Update values for detritus
  efficiencies[species %in% detritus] <- 
    exp(-1.670) * exp(0.164 * temp.kT) / (1 + exp(-1.670) * exp(0.164 * temp.kT))
  # Update values for plants
  efficiencies[species %in% plants]  <- 
    exp(0.179) * exp(0.164 * temp.kT) / (1 + exp(0.179) * exp(0.164 * temp.kT))
    
  # Check if some species forage on detritus and animals
  on.source <- 
    colSums(web_matrix[sources, ,drop = FALSE]) > 0
  on.anims <- 
    colSums(web_matrix[animals, ,drop = FALSE]) > 0
  if (sum(on.source & on.anims) > 0){
    cat("Warning: in web",site, web, "some species are feeding on detritus and animals: \n")
    cat(species[on.source & on.anims], '\n')
    if (any(colSums(web_matrix[non.basals,species[on.source & on.anims], drop = FALSE]) > 0.0)){
      cat("!!!!!  and some of the species are not feeding on fpom, cpom and detritus: \n")
      cat(species[on.source & on.anims & 
                    any(colSums(web_matrix[non.basals,species[on.source & on.anims], drop = FALSE]) > 0.0)], '\n')
    }
  }
  
  # Check if some species are feeding on the two detritus
  both.dets <- 
    colSums(web_matrix[detritus, ,drop = FALSE]) > 1
  if (sum(both.dets)>0){
    cat("Some species are eating on the two detritus: \n")
    cat(species[both.dets], '\n')}

  # Check if some species feed on plants and detritus
  if (any(colSums(web_matrix[plants, ,drop = FALSE]) > 0 & 
          colSums(web_matrix[detritus, ,drop = FALSE])>0)){
    cat("Some species are feeding on plant and detritus: \n")
    cat(species[colSums(web_matrix[plants, ,drop = FALSE]) > 0 & 
                  colSums(web_matrix[detritus, ,drop = FALSE])>0])
  }
  
  # Here adding a trycatch in case algorithm cannot compute the fluxes
  tryCatch(expr = {
    # Get fluxes
    mat.flux <-  
     fluxweb::fluxing(web_matrix, 
                      biomasses = biomasses, 
                      losses = met.rates, 
                      efficiencies = efficiencies)
    # Write as csv
   write.table(x = mat.flux, 
               sep = ",",
               file = paste0(here::here("data", "fluxes"), "/",
                             site, "_", web, ".txt"))},
   error = function(e) {
     print(paste("bromeliad", web, "in", site, "with problem"))
           })
}

# Prepare data for analysis ---------------------------------------------------------------
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
traits_CPI <-
  readr::read_csv(here::here("data",
                             "traits_CPI.csv"))

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
                              names = c("site", "bromeliad_id"))

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
  dplyr::left_join(traits_CPI %>% 
                     dplyr::select(dplyr::everything()) %>% 
                     dplyr::rename_with(~paste0("cons_",.x )) %>% 
                     dplyr::rename(consumer = cons_bwg_name),
                   by = "consumer") %>% 
  ## Add biomass data relative to consumers
  dplyr::left_join(readr::read_csv(here::here("data",
                                              "biomass_stoich.csv")) %>% 
                     dplyr::select(consumer = bwg_name, site = site_code ,bromeliad_id,
                                   tot_biomass_g, avg_bodymass_g, abundance, 
                                   cons_C = C, cons_N = N, cons_P = P) %>% 
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
  dplyr::left_join(traits_CPI %>% 
                     dplyr::select(bwg_name, functional_group, CPI_crisp) %>% 
                     dplyr::rename_with(~paste0("cons_",.x )) %>% 
                     dplyr::rename(consumer = cons_bwg_name),
                   by = "consumer") %>% 
  ## Add the composition of resource
  dplyr::left_join(traits_CPI %>% 
                     dplyr::select(bwg_name, functional_group, CPI_crisp) %>% 
                     dplyr::rename_with(~paste0("res_",.x )) %>% 
                     dplyr::rename(resource = res_bwg_name),
                   by = "resource") %>% 
  ## Add nutrient content for consumer
  dplyr::left_join(readr::read_csv(here::here("data",
                                              "biomass_stoich.csv")) %>% 
                     dplyr::select(consumer = bwg_name, site = site_code, bromeliad_id,
                                   cons_C = C, cons_N = N, cons_P = P) %>% 
                     dplyr::mutate(bromeliad_id = as.character(bromeliad_id)),
                   by = c("site", "bromeliad_id", "consumer")) %>% 
  ## Add nutrient content for resource
  dplyr::left_join(readr::read_csv(here::here("data",
                                              "biomass_stoich.csv")) %>% 
                     dplyr::select(resource = bwg_name, site = site_code, bromeliad_id,
                                   res_C = C, res_N = N, res_P = P) %>% 
                     dplyr::mutate(bromeliad_id = as.character(bromeliad_id)),
                   by = c("site", "bromeliad_id", "resource")) %>% 
  ## Get CN distance
  dplyr::mutate(CN_dist = sqrt((cons_C - res_C)^2 + (cons_N - res_N)^2)) %>% 
  ## Get CNP distance
  dplyr::mutate(CNP_dist = sqrt((cons_C - res_C)^2 + (cons_N - res_N)^2 + (cons_P - res_P)^2)) %>% 
  ## Replace distances of 0 with Nas
  dplyr::mutate(CN_dist = ifelse(CN_dist == 0, 
                                 NA, CN_dist),
                CNP_dist = ifelse(CNP_dist == 0, 
                                  NA, CNP_dist)) %>% 
  ## Rename flux_prop column
  dplyr::rename(flux_prop_pair = flux_prop,
                flux_raw_pair = flux_raw)



# Get Gini coefficient and total energy flux for each bromeliad, 
# based on all the fluxes (i.e. not pooled per consumer, excluding zooplankton)
fluxes_long_perbrom <- 
  fluxes_long_percons %>% 
  ## Remove zooplankton
  dplyr::filter(!stringr::str_detect(string = consumer,
                                     pattern = "plankton")) %>% 
  dplyr::select(site, bromeliad_id, flux_raw, cons_C, cons_N, cons_P) %>%
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
                     dplyr::select(site, bromeliad_id, flux_raw, cons_C, cons_N, cons_P) %>%
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

# Write data
readr::write_csv(fluxes_long_perbrom,
                 here::here("data",
                            "fluxes_long_perbrom.csv"))
readr::write_csv(fluxes_long_perpair,
                 here::here("data",
                            "fluxes_long_perpair_.csv"))
readr::write_csv(fluxes_long_percons,
                 here::here("data",
                            "fluxes_long_percons.csv"))



