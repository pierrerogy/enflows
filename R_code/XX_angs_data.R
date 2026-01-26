# Data from Ang's database
## https://datadryad.org/dataset/doi:10.5061/dryad.3tx95x6r2#readme

# Load libraries
library(tidyverse)
library(here)
library(readxl)
library(taxizedb)
source(here::here("R_code",
                  "functions.R"))


# Read names of sheets
sheet_list <-
  readxl::excel_sheets(here::here("data",
                                  "doi_10_5061_dryad_3tx95x6r2__v20250301",
                                  "datapaper_2025_02_11.xlsx"))

# Read all sheets into a list
og_data <-
  purrr::map(sheet_list, 
             ~readxl::read_excel(here::here("data",
                                            "doi_10_5061_dryad_3tx95x6r2__v20250301",
                                            "datapaper_2025_02_11.xlsx"), 
                                 sheet = .x)) %>%
  ## Add names
  purrr::set_names(sheet_list)

# Load data
dats <- 
  og_data$data %>% 
  ## Filter aquatic animals of interest
  dplyr::filter(habitat == "freshwater" & kingdom_init == "Animalia" &
                  phylum_init %in% c("Arthropoda", "Annelida",
                                     "Mollusca", "Rotifera",
                                     "Platyhelminthes", "Nematoda")) %>% 
  ### Remove that weird thing from the data
  dplyr::filter(data.source1 != "Rugenski 2013") 

# Let's check which species have changing taxonomy ------------------------
## Kingdom
dats %>% 
  dplyr::filter(kingdom_init != kingdom_revised)
### Nothing (PS I checked this before filtering in the last command)

## Phylum
dats %>% 
  dplyr::filter(phylum_init != phylum_revised)
### Nothing (PS I checked this before filtering in the last command)

## Class
dats %>% 
  dplyr::filter(class_init != class_revised) %>% 
  dplyr::select(data.source1,class_init, class_revised, 
                sp.morphosp_init, sp.morphosp_revised) %>% 
  unique()
### Two copepods

## Order
dats %>% 
  dplyr::filter(order_init != order_revised) %>% 
  dplyr::select(data.source1,order_init, order_revised,
                sp.morphosp_init, sp.morphosp_revised) %>% 
  unique()
### Annelids, Molluscs, Platyhelminthes, some typos w insects
### Some typos with insects: Leptonema was put Hemiptera, and Copelatus was put Syrphidae

## Family
dats %>% 
  dplyr::filter(family_init != family_revised) %>% 
  dplyr::select(data.source1,family_init, family_revised,
                sp.morphosp_init, sp.morphosp_revised) %>% 
  unique()
### Nothing major, small updates or subfamily becoming family
### Tipulidae to Limoniidae across several studies

## Species
dats %>% 
  dplyr::filter(sp.morphosp_init != sp.morphosp_revised) %>% 
  dplyr::select(data.source1,sp.morphosp_init, sp.morphosp_revised) %>% 
  unique() 
### Looks like a lot have punctuation added

## Put Tipulidae and Limoniidae together
dats <- 
  dats %>% 
  dplyr::mutate(family_revised == ifelse(family_revised %in% c("Tipulidae", "Limoniidae"),
                                        "Tipulidae_Limoniidae",
                                        family_revised))

# Need a genus column... --------------------------------------------------
# Make a new column
dats <- 
  dats %>% 
  ## First split column at space
  tidyr::separate(sp.morphosp_revised,
                  into = c("genus", "species"),
                  sep = " ",
                  remove = F) %>% 
  ## Some warnings, one double species name and two "non tanypodinae" but that's ok
  ## Remove species column
  dplyr::select(-species)

# Get vector of genera
genera <- 
  unique(c(dats$genus))

# Get ID for all
genera_id <- 
  taxizedb::name2taxid(genera,
                       out_type = "summary")

# Check which ones were removed
genera[which(genera %notin% genera_id$name)]
## Some genera were missed
## "Tentrepohlia", "Rhaphium", "Neocordulia", "Psephensus", "Dactylobaetis",
## "Megapodagrion", "Heterogrion" 
## One with typo: Cullicoides instead of Culicoides


# Add taxon rank and fix errors
genera_id <- 
  genera_id %>% 
  dplyr::mutate(rank = taxizedb::taxid2rank(genera_id$id)) %>% 
  dplyr::filter(rank == "genus") %>%
  dplyr::rename(genus = name) %>% 
  dplyr::select(-id, -rank) %>% 
  dplyr::pull() %>% 
  ## Add missing taxa
  c("Tentrepohlia", "Rhaphium", "Neocordulia", 
    "Psephensus", "Dactylobaetis", "Megapodagrion", 
    "Heterogrion", "Cullicoides")


# Fix genus column in data
dats <- 
  dats %>% 
  dplyr::mutate(genus = ifelse(genus %in% genera_id,
                               genus,
                               NA)) %>%
  ## Fix typo
  dplyr::mutate(genus = ifelse(genus == "Cullicoides",
                 "Culicoides",
                 genus)) %>% 
  ## Keep columns of interest
  dplyr::select(dplyr::contains("revised"), genus,
                C = taxa.c, N = taxa.n, P = taxa.p, ## % CNP
                bodymass_g = body.weight.gr) %>% 
  ## Give better names
  dplyr::rename_with(~ str_remove(., "_revised"), 
                     dplyr::everything()) %>% 
  dplyr::rename(species = sp.morphosp) %>% 
  ## And just because I like things organised
  dplyr::relocate(genus, 
                   .before = species) %>% 
  ## Last fix: subgenus Microculex has Culex for genus
  dplyr::mutate(genus = ifelse(stringr::str_detect(string = species,
                                                   pattern = "Microculex"),
                               "Culex",
                               genus))



# Get models best explaining CNP for different taxa ---------------------
# One more problem species in the data
dats <- 
  dats %>% 
  dplyr::filter(species != "Pelocoris sp.1")

# Split the data in a list with different taxonomic levels of interest - lm
model_pars <- 
  purrr::map(
    c("order", "family", "genus", "species"),
    ~ dats %>%
      ## Select only the level and size_mm
      dplyr::select(name = dplyr::all_of(.x), 
                    C, N, P, bodymass_g) %>% 
      ## Remove NAs, name being the taxonomic level
      dplyr::filter(!is.na(name) ) %>% 
      ## Add column with level name
      dplyr::mutate(level = .x)) %>%
  ## Now compute the models
  purrr::map(
    ~ .x %>%
      ## Group by taxon
      dplyr::group_by(name, level) %>%
      ## Use the safe model function to return NAs across when needed
      dplyr::summarise(dplyr::across(c(C, N, P),
                                     ~ list(safe_lm(cur_column(), 
                                                    cur_data(),
                                                    mode = "lm"))),
                       ## Keep group names!
                       .groups = "keep")) %>% 
  ## Stick all levels together
  purrr::list_c() %>% 
  ## Make longer
  tidyr::pivot_longer(cols = c(C, N, P),
                      names_to = "element",
                      values_to = "model") %>%
  ## Remove NAs
  dplyr::filter(!is.na(model)) %>% 
  ## Filter models with bad or near-perfect fit
  dplyr::filter(purrr::map_dbl(model,
                               ~ summary(.x)$r.squared) <= 0.95 &
                  purrr::map_dbl(model,
                                 ~ summary(.x)$r.squared) >= 0.3) %>%
  ## Filter models based on significance
  dplyr::filter(purrr::map_dbl(model,
                               ~ anova(.x)$`Pr(>F)`[1]) < 0.05)


# Split the data in a list with different taxonomic levels of interest - lm2
model_pars2 <- 
  purrr::map(
    c("order", "family", "genus", "species"),
    ~ dats %>%
      ## Select only the level and size_mm
      dplyr::select(name = dplyr::all_of(.x), 
                    C, N, P, bodymass_g) %>% 
      ## Remove NAs, name being the taxonomic level
      dplyr::filter(!is.na(name) ) %>% 
      ## Add column with level name
      dplyr::mutate(level = .x)) %>%
  ## Now compute the models
  purrr::map(
    ~ .x %>%
      ## Group by taxon
      dplyr::group_by(name, level) %>%
      ## Use the safe model function to return NAs across when needed
               dplyr::summarise(dplyr::across(c(C, N, P),
                                              ~ list(safe_lm(cur_column(), 
                                                             cur_data(),
                                                             mode = "lm2"))),
                                ## Keep group names!
                                .groups = "keep")) %>% 
  ## Stick all levels together
  purrr::list_c() %>% 
  ## Make longer
  tidyr::pivot_longer(cols = c(C, N, P),
                      names_to = "element",
                      values_to = "model") %>%
  ## Remove NAs
  dplyr::filter(!is.na(model)) %>% 
  ## Filter models with bad or near-perfect fit
  dplyr::filter(purrr::map_dbl(model, 
                               ~ .x$rsquare) <= 0.95 &
                  purrr::map_dbl(model, 
                                 ~ .x$rsquare) >= 0.3) %>% 
  ## Filter models based on significance
  dplyr::filter(purrr::map_dbl(model, 
                               ~ .x$P.param) < 0.05)

# Split the data in a list with different taxonomic levels of interest - brms
model_pars3 <- 
  purrr::map(
    c("order", "family", "genus", "species"),
    ~ dats %>%
      ## Select only the level and size_mm
      dplyr::select(name = dplyr::all_of(.x), 
                    C, N, P, bodymass_g) %>% 
      ## Remove NAs, name being the taxonomic level
      dplyr::filter(!is.na(name) ) %>% 
      ## Add column with level name
      dplyr::mutate(level = .x)) %>%
  ## Now compute the models
  purrr::map(
    ~ .x %>%
      ## Group by taxon
      dplyr::group_by(name, level) %>%
      ## Use the safe model function to return NAs across when needed
      dplyr::summarise(dplyr::across(c(C, N, P),
                                     ~ list(safe_lm(cur_column(), 
                                                    cur_data(),
                                                    mode = "brms"))),
                       ## Keep group names!
                       .groups = "keep")) %>% 
  ## Stick all levels together
  purrr::list_c() %>% 
  ## Make longer
  tidyr::pivot_longer(cols = c(C, N, P),
                      names_to = "element",
                      values_to = "model") %>%
  ## Remove NAs
  dplyr::filter(!is.na(model)) 



  ## Filter models with bad or near-perfect fit
  dplyr::filter(purrr::map_dbl(model, 
                               ~ .x$rsquare) <= 0.95 &
                  purrr::map_dbl(model, 
                                 ~ .x$rsquare) >= 0.3) %>% 
  ## Filter models based on significance
  dplyr::filter(purrr::map_dbl(model, 
                               ~ .x$P.param) < 0.05)
  

             
# See how many of each has how many elements
model_pars %>% 
  dplyr::group_by(level, name) %>% 
  dplyr::tally()

# See histogram of rsquared
model_pars$model %>% 
  purrr::map(
    ~summary(.x)$r.squared) %>% 
  purrr::list_c() %>% 
  hist()

# Get list of levels with relationships
saveRDS(object = model_pars,
        file = here::here("data",
                          "sig_stoich_model.RDS"))





# Olivier data for loss function ------------------------------------------

cpi <- 
  c(99, 28.8, 27, 10, 21, 27, 13.5, 22.5, 17, 54, 70)

pb <- 
  c(6.3, 11.8, 16.1, 10.7, 45.5, 72, 10.4, 35, 10, 10.6, 13.6)

da <- 
  c(0.28, 0.12, 0.035, 0.14, 0.018, 0.023, 0.03, 0.18, 0.083, 0.063, 0.71)


plot(pb~cpi)
lm1 <- lm(pb~cpi)
summary(lm1)
anova(lm1)

plot(pb~da)
lm2 <- lm(pb~da)
summary(lm2)
anova(lm2)

