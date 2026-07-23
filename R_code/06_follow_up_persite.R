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




# Figures S2 - S16 ----------------------

# Load the per-site models -------------------------------------------------
# Point to the folder that holds one RDS per site
persite_folder <- 
  here::here("data",
             "persite")

# Load every site's suite into a single named list (names = site codes)
persite_models <- 
  ## List the per-site RDS files
  list.files(path = persite_folder,
             pattern = "^models_.*[.]RDS$",
             full.names = TRUE) %>% 
  ## Name each path by its site code (drop the "models_" prefix and ".RDS" suffix)
  purrr::set_names(nm = basename(.) %>% 
                     stringr::str_remove(pattern = "^models_") %>% 
                     stringr::str_remove(pattern = "[.]RDS$")) %>% 
  ## Read each file into memory
  purrr::map(.f = readr::read_rds)

# Describe every model in the suite ----------------------------------------
# One row per model, in fit_site_suite() order: list position, data grain,
# x/y variables, basic_line() scaling rule, response offset (matching the fit,
# so points sit on the line), and output figure name.
plot_spec <- 
  tibble::tribble(
    ~position, ~grain,      ~x_var,          ~y_var,          ~scale,          ~y_offset, ~figname,
    # Hypothesis 1 - raw flow entering a consumer node ~ each PC axis
    1L,        "percons",   "Axis.1",        "flux_raw",      "x_scale",       0.001,     "figs2",
    2L,        "percons",   "Axis.2",        "flux_raw",      "x_scale",       0.001,     "figs3",
    3L,        "percons",   "Axis.3",        "flux_raw",      "x_scale",       0.001,     "figs4",
    4L,        "percons",   "Axis.4",        "flux_raw",      "x_scale",       0.001,     "figs5",
    # Hypothesis 2 - raw flow in a consumer-resource pair ~ functional distance
    5L,        "perpair",   "func_dist",     "flux_raw_pair", "x_scale",       0.001,     "figs6",
    6L,        "perpair",   "dist_axis1",    "flux_raw_pair", "x_scale",       0.001,     "figs7",
    7L,        "perpair",   "dist_axis2",    "flux_raw_pair", "x_scale",       0.001,     "figs8",
    8L,        "perpair",   "dist_axis3",    "flux_raw_pair", "x_scale",       0.001,     "figs9",
    9L,        "perpair",   "dist_axis4",    "flux_raw_pair", "x_scale",       0.001,     "figs10",
    # Hypothesis 3 - whole-web flow ~ functional volume and distance
    10L,       "perbrom",   "prop_hull",     "totalflow_Jy",  "none",          0,         "figs11",
    11L,       "perbrom",   "avg_func_dist", "totalflow_Jy",  "x_scale",       0,         "figs12",
    12L,       "perbrom",   "prop_hull",     "predflow_Jy",   "none",          0.01,      "figs13",
    13L,       "perbrom",   "avg_func_dist", "predflow_Jy",   "x_scale",       0.01,      "figs14",
    14L,       "perbrom",   "prop_hull",     "gini",          "y_original",    0,         "figs15",
    15L,       "perbrom",   "avg_func_dist", "gini",          "both_original", 0,         "figs16")

# Bundle the three data grains so each spec row can grab the right one
grain_data <- 
  list(percons = fluxes_long_percons,
       perpair = fluxes_long_perpair,
       perbrom = fluxes_long_perbrom)

# Draw and save one supplementary figure per model ------------------------
# Walk over the spec, one row (i.e. one model, hence one figure) at a time
purrr::pwalk(.l = plot_spec,
             .f = function(position, grain, x_var, y_var, scale, y_offset, figname){
               ## Sanity check - announce which figure is being built
               print(figname)
               
               ## Grab the data grain this model was fit on
               dats_grain <- 
                 grain_data[[grain]]
               
               ## Build one panel per site (imap hands us the site code as .y)
               siteplots <- 
                 purrr::imap(.x = persite_models,
                             .f = function(site_suite, site){
                               ### Pull this site's model for the focal relationship
                               model <- 
                                 site_suite[[position]]
                               
                               ### Subset raw data to this site, offset the response
                               ### (matches the fit, avoids log(0))
                               dats_site <-
                                 dats_grain %>% 
                                 dplyr::filter(site_code == site) %>% 
                                 dplyr::mutate(!!y_var := .data[[y_var]] + y_offset)
                               
                               ### Quadratic term label + bold flag (pd > 99.9, rope < 5)
                               quad_sig <-
                                 title_sig(model = model,
                                           index = 3)

                               ### Panel: site + N title, quadratic term annotation
                               basic_line(model = model,
                                          dats = dats_site,
                                          index = 1,
                                          x_var = x_var,
                                          y_var = y_var,
                                          scale = scale) +
                                 ggtitle(paste0(site,
                                                "; N = ",
                                                nrow(model$data))) +
                                 annotate(geom = "text",
                                          x = -Inf,
                                          y = Inf,
                                          hjust = 0,
                                          vjust = 2,
                                          label = quad_sig$label,
                                          fontface = quad_sig$fontface) +
                                 theme(axis.title.y = element_text(margin = margin(r = 0,  l = 10),
                                                                   size = rel(0.8)))})
               
               ## Arrange the site panels on a 6 x 5 grid (as in the original supplement)
               plot_grid <- 
                 cowplot::plot_grid(plotlist = siteplots,
                                    nrow = 6,
                                    ncol = 5)
               
               ## Save to the plots folder, one PDF per relationship
               ggplot2::ggsave(plot_grid,
                               height = 15,
                               width = 15,
                               filename = here::here("plots",
                                                     paste0(figname,
                                                            ".pdf")))})
