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
# And ultimately modified by  which is the most proximate source of the data here
## Srivastava, D. S. et al. 2023. Geographical variation in the trait‐based assembly 
## patterns of multitrophic invertebrate communities. Functional Ecology 37: 73–86.


# Libraries
library(tidyverse)
library(here)
library(hellometry)
source(here::here("R_code",
                  "functions.R"))

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

# PCA without functional groups relating to food  -------------------------------------------
# Get traits
fw_traits <- 
  ## Get fw data
  readr::read_csv(here::here("data", 
                             "traits_pca_ready.csv")) %>% 
  ## Use code from Céréghino et al. (2018) paper to remove species with more than 7 NAs in traits
  dplyr::filter(species_id %notin% c(4076, 4461, 4681, 5201, 5206, 5211, 5421, 6066, 6231, 
                                     6236, 6531, 6771, 6776, 6781, 6786, 6791, 6801, 6926,
                                     6936, 6941, 6946, 6951, 7251, 7316, 7886, 8026)) 


# Rank trait values
ranked_traits <- 
  fw_traits %>% 
  ## Keep only bwg name and traits that do not have to do with feeding mode
  dplyr::select(dplyr::matches("AS|BS|DM|RE|RF|LO|RM|MD|CP|BF",
                               ignore.case = F)) %>% 
  ## Rank each column according to values
  dplyr::mutate(across(.cols = everything(),
                       .fns = ~rank(.x,
                                    na.last = "keep"))) 

# Run PCA
RES_pca <-
  SYNCSA::pca(ranked_traits)

# Quick look at the variables
RES_pca$eigenvalues
head(RES_pca$individuals)
RES_pca$variables

# Plot in a grid
plot_grid <- 
  cowplot::plot_grid(plotlist = c(plz_plot_my_pca(PCA = RES_pca,
                                                  axes = c(1,2),
                                                  type = "points",
                                                  dats = fw_traits),
                                  plz_plot_my_pca(PCA = RES_pca,
                                                  axes = c(1,2),
                                                  type = "arrows",
                                                  dats = fw_traits),
                                  plz_plot_my_pca(PCA = RES_pca,
                                                  axes = c(1,3),
                                                  type = "points",
                                                  dats = fw_traits),
                                  plz_plot_my_pca(PCA = RES_pca,
                                                  axes = c(1,3),
                                                  type = "arrows",
                                                  dats = fw_traits),
                                  plz_plot_my_pca(PCA = RES_pca,
                                                  axes = c(1,4),
                                                  type = "points",
                                                  dats = fw_traits),
                                  plz_plot_my_pca(PCA = RES_pca,
                                                  axes = c(1,4),
                                                  type = "arrows",
                                                  dats = fw_traits)),
                     ncol = 2,
                     nrow = 3,
                     rel_widths = c (1.2, 1))



# Save plot
ggplot2::ggsave(here::here("plots",
                           "figs1.pdf"),
                plot_grid,
                width = 13,
                height = 13.5)

# Plot axes differently
cowplot::plot_grid(plotlist = c(plz_plot_my_pca(PCA = RES_pca,
                                                axes = c(2,3),
                                                type = "points",
                                                dats = fw_traits),
                                plz_plot_my_pca(PCA = RES_pca,
                                                axes = c(2,3),
                                                type = "arrows",
                                                dats = fw_traits),
                                plz_plot_my_pca(PCA = RES_pca,
                                                axes = c(2,4),
                                                type = "points",
                                                dats = fw_traits),
                                plz_plot_my_pca(PCA = RES_pca,
                                                axes = c(2,4),
                                                type = "arrows",
                                                dats = fw_traits),
                                plz_plot_my_pca(PCA = RES_pca,
                                                axes = c(3,4),
                                                type = "points",
                                                dats = fw_traits),
                                plz_plot_my_pca(PCA = RES_pca,
                                                axes = c(3,4),
                                                type = "arrows",
                                                dats = fw_traits)),
                   ncol = 2,
                   nrow = 3,
                   rel_widths = c (1.2, 1))




# Get PCA scores for all species
axes.raw <- 
  RES_pca$individuals

# Add new label
axes <- 
  axes.raw[,1:4]
## Double check data was read correcly
head(axes)
dim(axes)

# Scale for mean zero and standard deviation one
axes <- 
  scale(axes, 
        center = TRUE, 
        scale = TRUE)

# Now give the real name of all the species
axes <- 
  axes %>%
  ## Convert to tibble
  tibble::as_tibble() %>% 
  ## Get species code from pre-scaling data
  dplyr::mutate(species_id = as.character(fw_traits$species_id)) %>% 
  ## Join species names
  dplyr::left_join(species_names,
                   by = "species_id") %>% 
  ## remove species id
  dplyr::select(-species_id)

# Get distance between each two pairs of species --------------------------
## This is actually pretty easy
dist_values <- 
  ## First get overall distances
  dist(axes %>% 
         ## Make the species names rowwnames
         tibble::column_to_rownames(var = "bwg_name"), 
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
  dplyr::filter(resource != consumer) %>% 
  ## Add scores for consumer and resource
  ### Consumers
  dplyr::left_join(axes %>%
                     ## Make species name a column
                     dplyr::rename(consumer = bwg_name) %>% 
                     ## to not confuse scores for resources and consumers
                     dplyr::rename_with(~stringr::str_c("con_", .), 
                                        dplyr::contains("Axis")),
                   by = "consumer") %>%
  ### Resources
  dplyr::left_join(axes %>%
                     ## Make species name a column
                     dplyr::rename(resource = bwg_name) %>% 
                     ## to not confuse scores for resources and consumers
                     dplyr::rename_with(~stringr::str_c("res_", .), 
                                        dplyr::contains("Axis")),
                   by = "resource") %>%
  ## Get distances for single axes
  dplyr::mutate(dist_axis1 = abs(res_Axis.1 - con_Axis.1),
                dist_axis2 = abs(res_Axis.2 - con_Axis.2),
                dist_axis3 = abs(res_Axis.3 - con_Axis.3),
                dist_axis4 = abs(res_Axis.4 - con_Axis.4)) %>% 
  ## Remove scores for consumer and resources
  dplyr::select(consumer, resource, distance, starts_with("dist_axis"))


# Get the hull hypervolume fron the overall dataset -----------------------
# NB main difference here we do not remove outliers
# Compute observed convex hull with custom function
hull <- 
  geometry::convhulln(axes %>% 
                        ## Make the species names rowwnames
                        tibble::column_to_rownames(var = "bwg_name")
                        ,"FA")
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
                     dplyr::rename(consumer = bwg_name),
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
# Per consumer
## Combine
fluxes_long_percons <- 
  readr::read_csv(here::here("data",
                             "fluxes_long_percons.csv")) %>% 
  dplyr::left_join(axes %>% 
                     ## Rename bwg_name column to join
                     dplyr::rename(consumer = bwg_name),
                   by = c("consumer"))
## Save
readr::write_csv(fluxes_long_percons,
                 here::here("data",
                            "fluxes_long_percons.csv")) 

# Per pair
## Combine
fluxes_long_perpair <- 
  readr::read_csv(here::here("data",
                             "fluxes_long_perpair.csv")) %>% 
  dplyr::left_join(dist_values %>% 
                     dplyr::select(consumer, resource, func_dist = distance, 
                                   dplyr::contains("dist_axis")),
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
                     dplyr::summarise(avg_func_dist = mean(func_dist, na.rm = TRUE)),
                   by = c("site_code", "bromeliad_id"))
## Save
readr::write_csv(fluxes_long_perbrom,
                 here::here("data",
                            "fluxes_long_perbrom.csv")) 
