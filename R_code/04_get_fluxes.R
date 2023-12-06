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
                          biomass_kind = "both",
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
                     dplyr::select(-replicate, -visit_id) %>% 
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

# Add microorganism data
## Let's just randomly say that microzooplankton and mesozooplantkon are 5 and 10 times bacteria
## And bromeliad tissue is always 0.001
## FPOM is total detritus * 0.005 and CPOM is * 0.01
## IS TOTAL DETRITUS G OR MG?????????
basal_biomass <- 
  bromeliads_predicted %>% 
  dplyr::select(visit_id, bromeliad_id, site, total_detritus_mg = fpom_mg_dry, 
                algae = alg_g_estimate, microorganisms = bact_g_estimate) %>% 
  ## CHANGE HERE THE SCALE ONCE WE KNOW UNIT OF TTOAL DETRITUS
  dplyr::mutate(FPOM = total_detritus_mg * 0.005,
                CPOM = total_detritus_mg * 0.01,
                microzooplankton = microorganisms * 5,
                mesozooplankton = microorganisms * 5,
                bromeliad.tissue = 1) %>% 
  ## Remove empty rows
  dplyr::filter(!is.nan(algae)) %>% 
  ## Add site names
  dplyr::left_join(codes %>% 
                     dplyr::select(site_code, visit_id) %>% 
                     dplyr::mutate(visit_id = as.character(visit_id)),
                   by = "visit_id") %>% 
  ## Here couple problematic sites: cusuco just not in the data, and nouragues2006 wihtout code
  dplyr::mutate(site_code = ifelse(site == "Nouragues2006",
                                   "NOUR", site_code)) %>% 
  dplyr::filter(!is.na(site_code)) %>% 
  ## Remove superfluous columns
  dplyr::select(-visit_id, -site, -total_detritus_mg) %>% 
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
  dplyr::mutate(bromeliad_id = as.character(bromeliad_id)) %>% 
  dplyr::bind_rows(basal_biomass)

# Get food webs -----------------------------------------------------------
#  Register all site-level matrices 
site_matrices <- 
  load_matrices(folder = here::here("data",
                                    "probabilistic_matrices"),
                ending = "_NM.txt") %>% 
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
                   by = c("resource", "bromeliad_id", "site_code"))



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
    c("bromeliad.tissue")
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
