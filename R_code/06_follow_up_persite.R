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

site_matrices <-  site_matrices[18]

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
# Build the path to the output folder that will hold one RDS per site
persite_folder <-
  here::here("data",
             "persite")
# Check whether that folder already exists
if (!dir.exists(persite_folder)) {
  # Create it (and any missing parents) only if it is not there yet
  dir.create(persite_folder,
             recursive = TRUE)
}

# Map over sites, fitting the whole suite and saving each site as it completes
names(site_matrices) %>%
  # Keep site names so each fitted suite can be tied back to its site
  purrr::set_names() %>%
  # Walk (rather than map) since we save to disk for side effect, not in memory
  purrr::walk(function(site){
    # Fit the whole model suite for this one site
    models_site <-
      fit_site_suite(site = site,
                     percons = fluxes_long_percons,
                     perpair = fluxes_long_perpair,
                     perbrom = fluxes_long_perbrom)
    # Write this site's models to its own RDS file in the persite folder
    readr::write_rds(models_site,
                     file.path(persite_folder,
                               # Name the file after the site, e.g. models_XXX.RDS
                               paste0("models_",
                                      site,
                                      ".RDS")))
  })