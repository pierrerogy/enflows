# Get the functional coverage and distance for bromeliads

# Script adapted from script 2 in:
# Vanderlei J. Debastiani, Régis Céréghino, & Valério D. Pillar. (2018). 
# Code and results for Constraints on the functional trait space of aquatic 
# invertebrates in bromeliads (v1.0). Zenodo. https://doi.org/10.5281/zenodo.1200194
# Itself the data and code from:
# Céréghino, R. et al. 2018. Constraints on the functional trait space of aquatic 
# invertebrates in bromeliads. Funct Ecol 32: 2435–2447.
# Based on Diaz et al 2016 Nature:
# Díaz, S., Kattge, J., Cornelissen, J. H. C., Wright, I. J., Lavorel, S.,
# Dray, S., … Gorné, L. D. (2016). The global spectrum of plant form
# and function. Nature, 529, 167–171. https://doi.org/10.1038/nature16489
# ftp://pbil.univ-lyon1.fr/pub/datasets/dray/Diaz_Nature/

# Libraries
library(tidyverse)
library(here)

# Get list of species per bromeliad
species_list <- 
  readr::read_csv(here::here("data",
                            "fluxes_long_percons.csv")) %>% 
  dplyr::select(bromeliad_id, consumer)

# Get overall list of species names
species_names <- 
  readr::read_rds(here::here("data",
                             "scalewebs",
                             "fw_database.rds")) %>% 
  purrr::pluck("traits") %>% 
  dplyr::select(species_id, bwg_name) %>% 
  dplyr::distinct()

# Read in data (same file name than Zenodo repo)
axes.raw <- 
  read.table(here::here("data",
                        "RES_pca_individuals_0.7.7_ranks.txt"), 
             h = TRUE, row.names = 1)

# Add new label
axes <- 
  axes.raw[,1:4]
## Double check data was read correcly
head(axes)
dim(axes)

# Scale for mean zero and standard deviation one
axes <- 
  scale(axes, center = TRUE, scale = TRUE)

# Now give the real name of all the species
axes <- 
  axes %>%
  ## Convert to tibble
  tibble::as_tibble() %>% 
  ## Get species code from pre-scaling data
  dplyr::mutate(species_id = row.names(axes.raw)) %>% 
  ## Join species names
  dplyr::left_join(species_names,
                   by = "species_id") %>% 
  ## remove species id
  dplyr::select(-species_id) %>% 
  ## Make the species names rowwnames
  tibble::column_to_rownames(var = "bwg_name")
  

# Get distance between each two pairs of species --------------------------
## This is actually pretty easy
dist_values <- 
  ## Firs get distances
  dist(axes, 
       method = "euclidean") %>% 
  ## Convert to matrix
  as.matrix() %>% 
  ## Convert to long format
  as.data.frame() %>%
  ## Make row names resource column
  tibble::rownames_to_column(var = "resource") %>%
  ## Make data long format
  tidyr::pivot_longer(cols = where(is.numeric),
                      values_to = "distance",
                      names_to = "consumer") %>% 
  ## Remove cannibalism
  dplyr::filter(resource != consumer)


# Get the hull hypervolume fron the overall dataset -----------------------
# NB main difference here is that we do not remove outliers, whereas
# references did keep 95% of the data around the centroid

# Compute observed convex hull with custom function
hull <- 
  geometry::convhulln(axes,"FA")
## Extract the volume
overall_volume <- 
  hull$vol


# Get hypervolume for each bromeliad --------------------------------------

# List bromeliads food webs
# we will use basal from script 02
hypervolumes <- 
  species_list %>% 
  ## Remove basal species
  dplyr::filter(consumer %notin% basal) %>% 
  ## Add axes values
  dplyr::left_join(axes %>% 
                     ## Make species name a column
                     tibble::rownames_to_column(var = "consumer"),
                   by = "consumer") %>%
  ## Remove consumer column
  dplyr::select(-consumer) %>% 
  ## We have about 160 NAs because of species (mostly Coleoptera.20) with unknown traits, remove
  dplyr::filter(!is.na(Axis.1)) %>% 
  ## Convert to a list where each element is a bromeliad
  dplyr::group_by(bromeliad_id) %>%
  tidyr::nest(.key = "hull") %>% 
  ## Remove all bromeliads with less than 5 unique species (required for hull calculation)
  ### Get how many unique species
  dplyr::mutate(nspecies = purrr::map_int(hull,
                                          ~.x %>% 
                                            dplyr::distinct() %>% 
                                            nrow())) %>%
  ### Filter out and remove
  dplyr::filter(nspecies >= 5) %>% 
  dplyr::select(-nspecies) %>% 
  ## Compute hypervolume for each bromeliad
  dplyr::mutate(hull = purrr::map(hull,
                                  ~.x %>% 
                                    geometry::convhulln(.,
                                                        "FA"))) %>% 
  ## Only keep volume
  dplyr::mutate(hull = purrr::map(hull,
                                  ~.x$vol)) %>% 
  ## Flatten hull column
  dplyr::mutate(hull = unlist(hull)) %>% 
  ## Get proportion of overall hull
  dplyr::mutate(prop_hull = hull / overall_volume)
  


# Combine with data and save ----------------------------------------------
# Per pair
## Combine
fluxes_long_perpair <- 
  readr::read_csv(here::here("data",
                             "fluxes_long_perpair.csv")) %>% 
  dplyr::left_join(dist_values %>% 
                     dplyr::select(consumer, resource, func_dist = distance),
                   by = c("consumer", "resource"))
## Save
readr::write_csv(fluxes_long_perpair,
                 here::here("data",
                            "fluxes_long_perpair.csv")) 

# Per bromeliad
## Combine
fluxes_long_perbrom <- 
  readr::read_csv(here::here("data",
                             "fluxes_long_perbrom.csv")) %>% 
  dplyr::left_join(hypervolumes,
                   by = c("bromeliad_id")) %>% 
  ## Add average functional distance
  dplyr::left_join(fluxes_long_perpair %>% 
                     dplyr::select(func_dist, site_code, bromeliad_id) %>% 
                     dplyr::group_by(site_code, bromeliad_id) %>%
                     dplyr::summarise(avg_func_dist = mean(func_dist)),
                   by = c("site_code", "bromeliad_id"))
## Save
readr::write_csv(fluxes_long_perbrom,
                 here::here("data",
                            "fluxes_long_perbrom.csv")) 
