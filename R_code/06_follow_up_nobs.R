# Follow up analyses
# Load libraries 
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(bayestestR)
library(performance)

# Run the 4 MCMC chains of each model in parallel, one per allocated core
# (reads the cores given by SLURM, defaults to 4 for interactive runs)
options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                         unset = "4")))

# Load site matrices
site_matrices <- 
  load_matrices(folder = here::here("data",
                                    "matrices"),
                ending = "_PM.txt") %>% 
  make_binary_matrices(matrix_list = .,
                       cutoff = 0.02)

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

# Recompute trait PCA without BS traits -----------------------------------
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
  ## Keep only bwg name and traits that do not have to do with feeding mode, dropping BS
  dplyr::select(dplyr::matches("AS|DM|RE|RF|LO|RM|MD|CP|BF",
                               ignore.case = F)) %>%
  ## Rank each column according to values
  dplyr::mutate(across(.cols = everything(),
                       .fns = ~rank(.x,
                                    na.last = "keep")))

# Run PCA
RES_pca_nobs <-
  SYNCSA::pca(ranked_traits)

# # Visualise the PCA axes --------------------------------------------------
# # Quick look at the variables
# RES_pca_nobs$eigenvalues
# head(RES_pca_nobs$individuals)
# RES_pca_nobs$variables

# # Plot the first four axes (points and arrows) in a grid, as in script 04
# plot_grid_nobs <-
#   cowplot::plot_grid(plotlist = c(plz_plot_my_pca(PCA = RES_pca_nobs,
#                                                   axes = c(1,2),
#                                                   type = "points",
#                                                   dats = fw_traits),
#                                   plz_plot_my_pca(PCA = RES_pca_nobs,
#                                                   axes = c(1,2),
#                                                   type = "arrows",
#                                                   dats = fw_traits),
#                                   plz_plot_my_pca(PCA = RES_pca_nobs,
#                                                   axes = c(1,3),
#                                                   type = "points",
#                                                   dats = fw_traits),
#                                   plz_plot_my_pca(PCA = RES_pca_nobs,
#                                                   axes = c(1,3),
#                                                   type = "arrows",
#                                                   dats = fw_traits),
#                                   plz_plot_my_pca(PCA = RES_pca_nobs,
#                                                   axes = c(1,4),
#                                                   type = "points",
#                                                   dats = fw_traits),
#                                   plz_plot_my_pca(PCA = RES_pca_nobs,
#                                                   axes = c(1,4),
#                                                   type = "arrows",
#                                                   dats = fw_traits)),
#                      ncol = 2,
#                      nrow = 3,
#                      rel_widths = c (1.2, 1))

# # Display the grid
# plot_grid_nobs

# # Save the grid
# ggplot2::ggsave(here::here("plots",
#                            "figs2.pdf"),
#                 plot_grid_nobs,
#                 width = 13,
#                 height = 13.5)

# Prepare data with no-BS axes --------------------------------------------
# Basal species
basal <-
  c("algae","FPOM","CPOM","microorganisms",
    "microzooplankton","mesozooplankton", "bromeliad.tissue")

# Get overall list of species names
species_names <-
  readr::read_rds(here::here("data",
                             "scalewebs",
                             "fw_database.rds")) %>%
  purrr::pluck("traits") %>%
  dplyr::select(species_id, bwg_name) %>%
  dplyr::distinct()

# Get list of species per bromeliad
species_list <-
  fluxes_long_percons %>%
  dplyr::select(bromeliad_id, consumer)

# Get PCA scores for all species, keeping first 4 axes
axes <-
  RES_pca_nobs$individuals[,1:4] %>%
  ## Scale for mean zero and standard deviation one
  scale(center = TRUE,
        scale = TRUE) %>%
  ## Convert to tibble
  tibble::as_tibble() %>%
  ## Get species code from pre-scaling data
  dplyr::mutate(species_id = as.character(fw_traits$species_id)) %>%
  ## Join species names
  dplyr::left_join(species_names,
                   by = "species_id") %>%
  ## remove species id
  dplyr::select(-species_id)

# Get distance between each two pairs of species
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
                     dplyr::rename(consumer = bwg_name) %>%
                     dplyr::rename_with(~stringr::str_c("con_", .),
                                        dplyr::contains("Axis")),
                   by = "consumer") %>%
  ### Resources
  dplyr::left_join(axes %>%
                     dplyr::rename(resource = bwg_name) %>%
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

# Get the hull hypervolume from the overall dataset
overall_volume <-
  geometry::convhulln(axes %>%
                        tibble::column_to_rownames(var = "bwg_name"),
                      "FA")$vol

# Get hypervolume for each bromeliad
hypervolumes <-
  species_list %>%
  ## Remove basal species
  dplyr::filter(consumer %notin% basal) %>%
  ## Add axes values
  dplyr::left_join(axes %>%
                     dplyr::rename(consumer = bwg_name),
                   by = "consumer") %>%
  ## Remove consumer column
  dplyr::select(-consumer) %>%
  ## Remove species with unknown traits
  dplyr::filter(!is.na(Axis.1)) %>%
  ## Convert to a list where each element is a bromeliad
  dplyr::group_by(bromeliad_id) %>%
  tidyr::nest(.key = "hull") %>%
  ## Remove all bromeliads with less than 5 unique species (required for hull calculation)
  dplyr::mutate(nspecies = purrr::map_int(hull,
                                          ~.x %>%
                                            dplyr::distinct() %>%
                                            nrow())) %>%
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

# Combine no-BS axes with the flux data
## Per consumer
fluxes_long_percons_nobs <-
  fluxes_long_percons %>%
  ## Drop old axes
  dplyr::select(-dplyr::contains("Axis")) %>%
  dplyr::left_join(axes %>%
                     dplyr::rename(consumer = bwg_name),
                   by = "consumer")
## Per pair
fluxes_long_perpair_nobs <-
  fluxes_long_perpair %>%
  ## Drop old distances
  dplyr::select(-func_dist, -dplyr::contains("dist_axis")) %>%
  dplyr::left_join(dist_values %>%
                     dplyr::select(consumer, resource, func_dist = distance,
                                   dplyr::contains("dist_axis")),
                   by = c("consumer", "resource"))
## Per bromeliad
fluxes_long_perbrom_nobs <-
  fluxes_long_perbrom %>%
  ## Drop old hull and distance
  dplyr::select(-dplyr::any_of(c("hull", "prop_hull", "avg_func_dist"))) %>%
  dplyr::left_join(hypervolumes,
                   by = c("bromeliad_id")) %>%
  ## Add average functional distance
  dplyr::left_join(fluxes_long_perpair_nobs %>%
                     dplyr::select(func_dist, site_code, bromeliad_id) %>%
                     dplyr::group_by(site_code, bromeliad_id) %>%
                     dplyr::summarise(avg_func_dist = mean(func_dist, na.rm = TRUE)),
                   by = c("site_code", "bromeliad_id"))

# Hypothesis 1 - models at the scale of consumers -------------------------
# Non-linear model of raw energy flow entering consumer node ~ Axis.1
consflow_model_axis1_nobs <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.1,
                    scale = FALSE, center = TRUE) +
              I(scale(Axis.1, scale = FALSE, center = TRUE)^2) +
              (1 + scale(Axis.1,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons_nobs)
## Check assumptions
plot(consflow_model_axis1_nobs)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis1_nobs)
performance::r2_bayes(consflow_model_axis1_nobs)

# Non-linear model of raw energy flow entering consumer node ~ Axis.2
consflow_model_axis2_nobs <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.2,
                    scale = FALSE, center = TRUE) +
              I(scale(Axis.2, scale = FALSE, center = TRUE)^2) +
              (1 + scale(Axis.2,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons_nobs)
## Check assumptions
plot(consflow_model_axis2_nobs)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis2_nobs)
performance::r2_bayes(consflow_model_axis2_nobs)

# Non-linear model of raw energy flow entering consumer node ~ Axis.3
consflow_model_axis3_nobs <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.3,
                    scale = FALSE, center = TRUE) +
              I(scale(Axis.3, scale = FALSE, center = TRUE)^2) +
              (1 + scale(Axis.3,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons_nobs)
## Check assumptions
plot(consflow_model_axis3_nobs)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis3_nobs)
performance::r2_bayes(consflow_model_axis3_nobs)

# Non-linear model of raw energy flow entering consumer node ~ Axis.4
consflow_model_axis4_nobs <-
  brms::brm(log(flux_raw + 0.001) ~
              scale(Axis.4,
                    scale = FALSE, center = TRUE) +
              I(scale(Axis.4, scale = FALSE, center = TRUE)^2) +
              (1 + scale(Axis.4,
                         scale = FALSE, center = TRUE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15),
            data = fluxes_long_percons_nobs)
## Check assumptions
plot(consflow_model_axis4_nobs)
## Look at the effect
bayestestR::describe_posterior(consflow_model_axis4_nobs)
performance::r2_bayes(consflow_model_axis4_nobs)

# Save models
readr::write_rds(x = consflow_model_axis1_nobs,
                 file = here::here("data",
                                   "consflow_model_axis1_nobs.RDS"))
readr::write_rds(x = consflow_model_axis2_nobs,
                 file = here::here("data",
                                   "consflow_model_axis2_nobs.RDS"))
readr::write_rds(x = consflow_model_axis3_nobs,
                 file = here::here("data",
                                   "consflow_model_axis3_nobs.RDS"))
readr::write_rds(x = consflow_model_axis4_nobs,
                 file = here::here("data",
                                   "consflow_model_axis4_nobs.RDS"))

# Hypothesis 2 - more functional dissimilarity in a pair means less energy flow ----
# Overall functional distance
pairflow_model_overall_nobs <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(func_dist,
                    center = TRUE, scale = FALSE) +
              I(scale(func_dist, center = TRUE, scale = FALSE)^2) +
              (1 + scale(func_dist,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair_nobs)
## Check assumptions
plot(pairflow_model_overall_nobs)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_overall_nobs)
performance::r2_bayes(pairflow_model_overall_nobs)

# Distance for Axis 1
pairflow_model_axis1_nobs <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis1,
                    center = TRUE, scale = FALSE) +
              I(scale(dist_axis1, center = TRUE, scale = FALSE)^2) +
              (1 + scale(dist_axis1,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair_nobs)
## Check assumptions
plot(pairflow_model_axis1_nobs)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis1_nobs)
performance::r2_bayes(pairflow_model_axis1_nobs)

# Distance for Axis 2
pairflow_model_axis2_nobs <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis2,
                    center = TRUE, scale = FALSE) +
              I(scale(dist_axis2, center = TRUE, scale = FALSE)^2) +
              (1 + scale(dist_axis2,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair_nobs)
## Check assumptions
plot(pairflow_model_axis2_nobs)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis2_nobs)
performance::r2_bayes(pairflow_model_axis2_nobs)

# Distance for Axis 3
pairflow_model_axis3_nobs <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis3,
                    center = TRUE, scale = FALSE) +
              I(scale(dist_axis3, center = TRUE, scale = FALSE)^2) +
              (1 + scale(dist_axis3,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair_nobs)
## Check assumptions
plot(pairflow_model_axis3_nobs)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis3_nobs)
performance::r2_bayes(pairflow_model_axis3_nobs)

# Distance for Axis 4
pairflow_model_axis4_nobs <-
  brms::brm(log(flux_raw_pair + 0.001) ~
              scale(dist_axis4,
                    center = TRUE, scale = FALSE) +
              I(scale(dist_axis4, center = TRUE, scale = FALSE)^2) +
              (1 + scale(dist_axis4,
                         center = TRUE, scale = FALSE)|site_code) +
              (1|site_code:bromeliad_id),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perpair_nobs)
## Check assumptions
plot(pairflow_model_axis4_nobs)
## Look at the effect
bayestestR::describe_posterior(pairflow_model_axis4_nobs)
performance::r2_bayes(pairflow_model_axis4_nobs)

# Save models
readr::write_rds(x = pairflow_model_overall_nobs,
                 file = here::here("data",
                                   "pairflow_model_overall_nobs.RDS"))
readr::write_rds(x = pairflow_model_axis1_nobs,
                 file = here::here("data",
                                   "pairflow_model_axis1_nobs.RDS"))
readr::write_rds(x = pairflow_model_axis2_nobs,
                 file = here::here("data",
                                   "pairflow_model_axis2_nobs.RDS"))
readr::write_rds(x = pairflow_model_axis3_nobs,
                 file = here::here("data",
                                   "pairflow_model_axis3_nobs.RDS"))
readr::write_rds(x = pairflow_model_axis4_nobs,
                 file = here::here("data",
                                   "pairflow_model_axis4_nobs.RDS"))

# Hypothesis 3 - more and less imbalanced flow if less functional distances in whole food web ----
# Total flow
## Volume
totalflow_model_vol_nobs <-
  brms::brm(log(totalflow_Jy) ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
              I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom_nobs)
## Check assumptions
plot(totalflow_model_vol_nobs)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_vol_nobs)
performance::r2_bayes(totalflow_model_vol_nobs)
## Distance
totalflow_model_dist_nobs <-
  brms::brm(log(totalflow_Jy) ~
              scale(avg_func_dist,
                    center = TRUE, scale = FALSE) +
              I(scale(avg_func_dist, center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom_nobs)
## Check assumptions
plot(totalflow_model_dist_nobs)
## Look at the effect
bayestestR::describe_posterior(totalflow_model_dist_nobs)
performance::r2_bayes(totalflow_model_dist_nobs)

# Flow to predators
## Volume
predflow_model_vol_nobs <-
  brms::brm(log(predflow_Jy + 0.01) ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
              I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom_nobs)
## Check assumptions
plot(predflow_model_vol_nobs)
## Look at the effect
bayestestR::describe_posterior(predflow_model_vol_nobs)
performance::r2_bayes(predflow_model_vol_nobs)
## Distance
predflow_model_dist_nobs <-
  brms::brm(log(predflow_Jy + 0.01) ~
              scale(avg_func_dist,
                    center = TRUE, scale = FALSE) +
              I(scale(avg_func_dist, center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom_nobs)
## Check assumptions
plot(predflow_model_dist_nobs)
## Look at the effect
bayestestR::describe_posterior(predflow_model_dist_nobs)
performance::r2_bayes(predflow_model_dist_nobs)

# Gini index
## Volume
gini_model_vol_nobs <-

  brms::brm(gini ~
              scale(log(prop_hull),
                    center = TRUE, scale = FALSE) +
              I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom_nobs)
## Check assumptions
plot(gini_model_vol_nobs)
## Look at the effect
bayestestR::describe_posterior(gini_model_vol_nobs)
performance::r2_bayes(gini_model_vol_nobs)
## Distance
gini_model_dist_nobs <-
  brms::brm(gini ~
              scale(log(avg_func_dist),
                    center = TRUE, scale = FALSE) +
              I(scale(log(avg_func_dist), center = TRUE, scale = FALSE)^2) +
              (1 + scale(log(prop_hull),
                         center = TRUE, scale = FALSE) |site_code),
            iter = 20000,
            family = gaussian(),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = fluxes_long_perbrom_nobs)
## Check assumptions
plot(gini_model_dist_nobs)
## Look at the effect
bayestestR::describe_posterior(gini_model_dist_nobs)
performance::r2_bayes(gini_model_dist_nobs)

# Save models
readr::write_rds(x = totalflow_model_vol_nobs,
                 file = here::here("data",
                                   "totalflow_model_vol_nobs.RDS"))
readr::write_rds(x = predflow_model_vol_nobs,
                 file = here::here("data",
                                   "predflow_model_vol_nobs.RDS"))
readr::write_rds(x = gini_model_vol_nobs,
                 file = here::here("data",
                                   "gini_model_vol_nobs.RDS"))
readr::write_rds(x = totalflow_model_dist_nobs,
                 file = here::here("data",
                                   "totalflow_model_dist_nobs.RDS"))
readr::write_rds(x = predflow_model_dist_nobs,
                 file = here::here("data",
                                   "predflow_model_dist_nobs.RDS"))
readr::write_rds(x = gini_model_dist_nobs,
                 file = here::here("data",
                                   "gini_model_dist_nobs.RDS"))
