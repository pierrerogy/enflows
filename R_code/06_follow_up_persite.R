# Follow up analyses - same suite of models as the main analysis, fitted site-by-site
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

# Load site matrices (used here only to get the list of sites)
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

# Fit the suite for every site --------------------------------------------
# Map over sites, returning a list of lists (one named list of models per site)
models_persite <-
  names(site_matrices) %>%
  ## Keep site names so the output list is named
  purrr::set_names() %>%
  ## Fit the whole suite for each site
  purrr::map(~ fit_site_suite(site = .x,
                              percons = fluxes_long_percons,
                              perpair = fluxes_long_perpair,
                              perbrom = fluxes_long_perbrom))

# Save the list of lists
readr::write_rds(models_persite,
                 here::here("data",
                            "models_persite.RDS"))