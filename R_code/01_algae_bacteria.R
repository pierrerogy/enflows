# Prepare data
# Author: Pierre Rogy

# Libraries
library(here)
library(tidyverse)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(bayestestR)

# Get chelsa data (txt file obtained from the chelsa website)
get_chelsa(txt_file = here::here("data",
                                 "envidatS3paths.txt"),
           output_folder = here::here("chelsa"))

# Prepare drought experiment data -------------------------------------------------
# French Guyana
dats1 <- 
  readr::read_csv(here::here("data",
                           "closed area.csv"),
                  locale=locale(encoding="latin1")) %>% 
  ## Select and rename relevant columns
  dplyr::select(brom_sp = `Bromeliad sp`,
                brom_ID = `Brome ID`,
                trt.name,
                chl_a_ugL = `Chl a µg/L`,
                alg_mL = `Algae/mL`,
                bact_mL = `total bact/ml`,
                diam1 = `Diameter1 (cm)`,
                diam2 = `Diameter2 (cm)`,
                n_green_leaves = `# green Leaves`,
                actual_water_mL = `Actual Water Volume (ml)`,
                fpom_mg_dry = `FPOM (mg dry weight)`,
                incident_radiation_prcnt = `% incident radiation`,
                water_holding_cap_mL = maxvol) %>% 
  ## Make new columns
  dplyr::mutate(mean_diameter = (diam1 + diam2)/ 2,
                site = "nouragues",
                lat = 4.083,
                lon = -52.6830,
                brom_ID = as.character(brom_ID)) %>% 
  ## Remove old columns
  dplyr::select(-diam1, - diam2)

## Las Gamas and Macae
dats2 <- 
  readr::read_csv(here::here("data",
                             "argbra.csv")) %>% 
  ## Make new columns and replace negative values with 0
  dplyr::mutate(lat = ifelse(site == "las_gamas", 
                             -29.433, -22.2075),
                lon = ifelse(site == "las_gamas", 
                             -60.467, -41.4847),
                chl_a_ugL = ifelse(chl_a_ugL < 0, 
                                   0, chl_a_ugL)) 
# Combine data
dats_drought <- 
  dats2 %>% 
  dplyr::bind_rows(dats1)

# Add environmental data
dats_drought <- 
  add_env(dats = dats_drought,
          raster_folder = here::here("chelsa"))


# French Guyana data ------------------------------------------------------
## Régis
# voilà ce que l’on a sur des broméliacées non manipulées, c’est à dire de
# s données d’observation hors expérimentations. Pour ce qui est de l’ouverture de
# la canopée, nous mesurons le % de radiation incidente avec des photos 
# hémisphériques puis ce pourcentage est moyenné sur l’année à l’aide des co
# ordonnées GPS de chaque broméliacée (on calcule ça avec le logiciel Gap Light
# Analyser). Globalement, ce sont plutôt des sites ouverts ) semi-ouverts.
dats_fg <- 
  readr::read_csv(here::here("data",
                             "Data_Algae_Bact_FG.csv")) %>% 
  dplyr::rename(year = Year,
                site = Site,
                brom_ID = `plant ID`,
                brom_sp = `Plant species`,
                actual_water_mL = `Actual Water volume (ml)`,
                n_leaves = `N green leaves`,
                incident_radiation_prcnt = `Incident radiation (%)`,
                fpom_mg_dry = `FPOM <150 (mg dry weight)`,
                bact_mL = `Bacteria_number/ml`,
                alg_mL =`Algae number/ml`) %>% 
  ## Make new columns and replace negative values with 0
  dplyr::mutate(lat = ifelse(site == "Kaw", 
                             4.5146, ifelse(site == "Rochambeau",
                                            4.48849, ifelse(site == "Petit-Saut",
                                                            5.0629, ifelse(site == "Nouragues",
                                                                           4.083, ifelse(site == "Virginie",
                                                                                         4.18449, ifelse(site == "Sinnamary",
                                                                                                         5.3746, 5.24327)))))),
                lon = ifelse(site == "Kaw", 
                             -52.0665, ifelse(site == "Rochambeau",
                                               -52.21573, ifelse(site == "Petit-Saut",
                                                                 -53.0475, ifelse(site == "Nouragues",
                                                                                  -52.683, ifelse(site == "Virginie",
                                                                                                  -52.13543, ifelse(site == "Sinnamary",
                                                                                                                    -52.9546, -53.39197)))))),
                site = parse_character(site))
# Add environmental data
dats_fg <- 
  add_env(dats = dats_fg,
          raster_folder = here::here("chelsa"))

# Last tweaks on the data -------------------------------------------------
## Drought
equ_dats_drought <- 
  dats_drought %>% 
  ## Composite openness measurement
  dplyr::mutate(openness = mean_diameter/n_green_leaves) %>% 
  ## Log
  dplyr::mutate(across(c(chl_a_ugL, bact_mL,
                         actual_water_mL, fpom_mg_dry,openness,
                         incident_radiation_prcnt, npp,turbidity_final,
                         water_holding_cap_mL, bio1, bio15, ai),
                       ~as.numeric(log(.+0.001)),
                       .names = "{.col}_log")) %>% 
  ## Scale
  dplyr::mutate(across(c(chl_a_ugL, bact_mL,
                         actual_water_mL, fpom_mg_dry,openness,
                         incident_radiation_prcnt, npp, turbidity_final,
                         water_holding_cap_mL, bio1, bio15, ai),
                       ~as.numeric(scale(log(.+0.001))),
                       .names = "{.col}_log_scale")) %>% 
  ## convert incident radiation to open/close
  dplyr::mutate(canopy_binary = ifelse(incident_radiation_prcnt < 50,
                                       0, 1))

## FG  
equ_dats_fg <- 
  dats_fg %>% 
  ## Just rename one site
  dplyr::mutate(site = ifelse(site == "Angoul\xe8me",
                              "Angouleme", site)) %>% 
  ## Log
  dplyr::mutate(across(c(alg_mL, bact_mL,
                         actual_water_mL, fpom_mg_dry,
                         incident_radiation_prcnt,
                         bio1, bio15, ai, npp),
                       ~as.numeric(log(.+0.001)),
                       .names = "{.col}_log")) %>% 
  ## Scale
  dplyr::mutate(across(c(alg_mL, bact_mL,
                         actual_water_mL, fpom_mg_dry,
                         incident_radiation_prcnt, npp,
                         bio1, bio15, ai),
                       ~as.numeric(scale(log(.+0.001))),
                       .names = "{.col}_log_scale"))  %>% 
  ## convert incident radiation to open/close
  dplyr::mutate(canopy_binary = ifelse(incident_radiation_prcnt < 50,
                                       0, 1))

# Load fw data
fw_database <- 
  readr::read_rds(here::here("data", "other",
                             "fw_database.rds"))


# Algae - drought expt -------------------------------------------------------------------
# Make all possible combinations of variables
vars <- 
  ## Generate all possible combinations of variables
  tidyr::crossing(
    a = c(NA, "bact_mL_log"),
    b = c(NA, "actual_water_mL_log"),
    c = c(NA, "canopy_binary"),
    d = c(NA, "fpom_mg_dry_log"),
    e = c(NA, "bio1_log"),
    f = c(NA, "bio15_log"),
    g = c(NA, "ai_log")) %>% 
  ## Combine with a + to get model specs
  tidyr::unite(a, b, c, d, e, f, g,
               sep = " + ") %>% 
  ## Remove the NAS
  dplyr::mutate(a = stringr::str_remove_all(string = a,
                                            pattern = " [+] NA|NA [+] ")) %>% 
  ## Add the dependent variable and random effect
  dplyr::mutate(a = paste("chl_a_ugL_log ~", a, "+ (1|site)")) %>% 
  ## Remove the last NAs
  dplyr::filter(!stringr::str_detect(string = a,
                                     pattern = "~ NA| ~ [+] (1|site)")) %>% 
  ## Pull to make vector
  dplyr::pull() %>% 
  ## Keep unique combinations
  unique()

# Make list to store all models
alg_mods <- 
  list()

# Make loop that will make all models and save them in the list
for(i in 1:length(vars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(vars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = equ_dats_drought, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      alg_mods[[length(alg_mods) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(vars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = equ_dats_drought, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      alg_mods[[length(alg_mods) + 1]] <- 
        model
  }
  
}

# Get loo
alg_loo <- 
  alg_mods %>% 
  purrr::map(~loo(.x))

# Add names
names(alg_loo) <- 
  alg_mods %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()


## extract n of different points
npoints <- 
  alg_loo %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list
  assign(
    paste0("alg_loo_", i),
    alg_loo %>% 
      purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
      brms::loo_compare() %>% 
      get_eq_models()
  )
}

# Get most parsimonious model in each category
alg_loo_1 # (fpom_mg_dry_log)
alg_loo_2 # (actual_water_mL_log) | (ai_log) | (bio1_log))   


# Compute these models one by one and examine
## (fpom_mg_dry_log)
alg_mod1 <- 
  brms::brm(formula = chl_a_ugL_log ~ 
              fpom_mg_dry_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            data = equ_dats_drought, 
            family = gaussian())
## (actual_water_mL_log)
alg_mod2a <- 
  update(alg_mod1,
         formula = chl_a_ugL_log ~
           actual_water_mL_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_drought, 
         family = gaussian())
## (ai_log)
alg_mod2b <- 
  update(alg_mod1,
         formula = chl_a_ugL_log ~
           ai_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_drought, 
         family = gaussian())
## (bio1_log)
alg_mod2c <- 
  update(alg_mod1,
         formula = chl_a_ugL_log ~
           bio1_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_drought, 
         family = gaussian())

# Compare the effect of independent variables
bayestestR::describe_posterior(alg_mod1) # fpom sig
bayestestR::describe_posterior(alg_mod2a) # actual water sig
bayestestR::describe_posterior(alg_mod2b) # ai sig
bayestestR::describe_posterior(alg_mod2c) # bio1 not sig

# Compare the rsquared
brms::bayes_R2(alg_mod1) # 0.43
brms::bayes_R2(alg_mod2a) # 0.43
brms::bayes_R2(alg_mod2b) # 0.41
brms::bayes_R2(alg_mod2c) # 0.41

# Algae - FG -------------------------------------------------------------------
# Make all possible combinations of variables
vars <- 
  ## Generate all possible combinations of variables
  tidyr::crossing(
    a = c(NA, "bact_mL_log"),
    b = c(NA, "actual_water_mL_log"),
    c = c(NA, "canopy_binary"),
    d = c(NA, "fpom_mg_dry_log"),
    e = c(NA, "bio1_log"),
    f = c(NA, "bio15_log"),
    g = c(NA, "ai_log")) %>% 
  ## Combine with a + to get model specs
  tidyr::unite(a, b, c, d, e, f, g,
               sep = " + ") %>% 
  ## Remove the NAS
  dplyr::mutate(a = stringr::str_remove_all(string = a,
                                            pattern = " [+] NA|NA [+] ")) %>% 
  ## Add the dependent variable
  dplyr::mutate(a = paste("alg_mL_log ~", a, "+ (1|site)")) %>% 
  ## Remove the last NAs
  dplyr::filter(!stringr::str_detect(string = a,
                                     pattern = "~ NA| ~ [+] (1|site)")) %>% 
  ## Pull to make vector
  dplyr::pull() %>% 
  ## Keep unique combinations
  unique()

# Make list to store all models
alg_mods2 <- 
  list()

# Make loop that will make all models and save them in the list
for(i in 1:length(vars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(vars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = equ_dats_fg, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      alg_mods2[[length(alg_mods2) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(vars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = equ_dats_fg, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      alg_mods2[[length(alg_mods2) + 1]] <- 
        model
  }
  
}

# Get loo
alg_loo2 <- 
  alg_mods2 %>% 
  purrr::map(~loo(.x))

# Add names
names(alg_loo2) <- 
  alg_mods2 %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()


## extract n of different points
npoints <- 
  alg_loo2 %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list
  assign(
    paste0("alg_loo_", i),
    alg_loo2 %>% 
      purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
      brms::loo_compare() %>% 
      get_eq_models()
  )
}

# Get most parsimonious model in each category
alg_loo_1 # (fpom_mg_dry_log) 
alg_loo_2 #  actual_water_mL_log | bio1_log | bio15_log | ai_log | canopy_binary

# Compute these models one by one and examine
## (fpom_mg_dry_log)
alg_mod3 <- 
  brms::brm(formula = alg_mL_log ~ 
              fpom_mg_dry_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_fg, 
            family = gaussian())


## (actual_water_mL_log)
alg_mod4a <- 
  update(alg_mod3a,
         formula = alg_mL_log ~ 
           actual_water_mL_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
## (bio1_log)
alg_mod4b <- 
  update(alg_mod3a,
         formula = alg_mL_log ~
           bio1_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
## (bio15_log)
alg_mod4c <- 
  update(alg_mod3a,
         formula = alg_mL_log ~
           bio15_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
## (ai_log)
alg_mod4d <- 
  update(alg_mod3a,
         formula = alg_mL_log ~
           ai_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
## (canopy_binary)
alg_mod4e <- 
  update(alg_mod3a,
         formula = alg_mL_log ~
           canopy_binary + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())

# Compare the effect of independent variables
bayestestR::describe_posterior(alg_mod3) # fpom not sig
bayestestR::describe_posterior(alg_mod4a) # actual water not sig
bayestestR::describe_posterior(alg_mod4b) # bio1 not sig
bayestestR::describe_posterior(alg_mod4c) # bio15 not sig
bayestestR::describe_posterior(alg_mod4d) # ai not sig
bayestestR::describe_posterior(alg_mod4e) # canopy binary not sig
# Compare the rsquared
brms::bayes_R2(alg_mod3) # 0.07
brms::bayes_R2(alg_mod4a) # 0.20
brms::bayes_R2(alg_mod4b) # 0.19
brms::bayes_R2(alg_mod4c) # 0.18
brms::bayes_R2(alg_mod4d) # 0.19
brms::bayes_R2(alg_mod4e) # 0.18

# Bact, both datasets together  ---------------------------------------------
# Microorganism measurment in FG is flux cytometry, drought expt method is similar
equ_dats_both <- 
  equ_dats_drought %>% 
  ## will have to rescale everything here
  dplyr::select(-dplyr::contains("log")) %>% 
  dplyr::bind_rows(equ_dats_fg  %>% 
                     ## will have to rescale everything here
                     dplyr::select(-dplyr::contains("log"))) %>% 
  ## log and scale
  ## Log
  dplyr::mutate(across(c(chl_a_ugL, bact_mL,alg_mL, 
                         actual_water_mL, fpom_mg_dry,openness,
                         incident_radiation_prcnt, npp, turbidity_final,
                         water_holding_cap_mL, bio1, bio15, ai),
                       ~as.numeric(log(.+0.001)),
                       .names = "{.col}_log")) %>% 
  ## Scale
  dplyr::mutate(across(c(chl_a_ugL, bact_mL, alg_mL,
                         actual_water_mL, fpom_mg_dry,openness,
                         incident_radiation_prcnt, npp, turbidity_final,
                         water_holding_cap_mL, bio1, bio15, ai),
                       ~as.numeric(scale(log(.+0.001))),
                       .names = "{.col}_log_scale"))

# Get list of vars
vars <- 
  ## Generate all possible combinations of variables
  tidyr::crossing(
    a = c(NA, "actual_water_mL_log"),
    b = c(NA, "fpom_mg_dry_log"),
    c = c(NA, "canopy_binary"),
    d = c(NA, "bio1_log"),
    e = c(NA, "bio15_log"),
    f = c(NA, "ai_log"),
    g = c(NA, "alg_mL_log"),
    h = c(NA, "chl_a_ugL_log")) %>% 
  ## Combine with a + to get model specs
  tidyr::unite(a, b, c, d, e, f, g, h, 
               sep = " + ") %>% 
  ## Remove the NAS or the two algae measuremnts when together
  dplyr::mutate(a = stringr::str_remove_all(string = a,
                                            pattern = " [+] NA|NA [+] | [+] alg_mL_log [+] chl_a_ugL_log")) %>% 
  ## Add the dependent variable
  dplyr::mutate(a = paste("bact_mL_log ~", a, "+ (1|site)")) %>% 
  ## Remove the last NAs
  dplyr::filter(!stringr::str_detect(string = a,
                                     pattern = "~ NA|~ alg_mL_log [+] chl_a_ugL_log| ~ [+] (1|site)")) %>% 
  ## Pull to make vector
  dplyr::pull() %>% 
  ## Keep unique combinations
  unique()

# Make list to store all models
bact_mods <- 
  list()

# Make loop that will make all models and save them in the list
for(i in 1:length(vars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(vars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = equ_dats_both, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      bact_mods[[length(bact_mods) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(vars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = equ_dats_both, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      bact_mods[[length(bact_mods) + 1]] <- 
        model
  }
}  

# Get loo
bact_loo <- 
  bact_mods %>% 
  purrr::map(~loo(.x))

# Add names
names(bact_loo) <- 
  bact_mods %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()

## extract n of different points
npoints <- 
  bact_loo %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list
  assign(
    paste0("bact_loo_", i),
    bact_loo %>% 
      purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
      brms::loo_compare() %>% 
      get_eq_models()
  )
}

# Get most parsimonious model in each category
bact_loo_1 # fpom_mg_dry_log +  alg_mL_log
bact_loo_2 # (fpom_mg_dry_log)  
bact_loo_3 # chl_a_ugL_log
bact_loo_4 # alg_mL_log
bact_loo_5 # (ai_log) | (bio1_log) | (bio15_log)

# Compute model
# LOO 1
bact_mods1 <- 
  brms::brm(formula = bact_mL_log ~ fpom_mg_dry_log + alg_mL_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())

# LOO 2
bact_mods2 <- 
  brms::brm(formula = bact_mL_log ~ fpom_mg_dry_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())

# LOO 3
bact_mods3 <- 
  brms::brm(formula = bact_mL_log ~ chl_a_ugL_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())

# LOO 4
bact_mods4 <- 
  brms::brm(formula = bact_mL_log ~ alg_mL_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())

# LOO 5
bact_mods5a <- 
  brms::brm(formula = bact_mL_log ~ bio15_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())
bact_mods5b <- 
  brms::brm(formula = bact_mL_log ~  bio1_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())
bact_mods5c <- 
  brms::brm(formula = bact_mL_log ~ ai_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())

# Filter models based on characteristics
# LOO 5 out because only site level variables
# Compare the effect of independent variables
bayestestR::describe_posterior(bact_mods1) # fpom_log  alg_log not sig
bayestestR::describe_posterior(bact_mods2) # fpom sig
bayestestR::describe_posterior(bact_mods3) # chloro not sig
bayestestR::describe_posterior(bact_mods4) # alg  sig
bayestestR::describe_posterior(bact_mods5a) # bio15 not sig
bayestestR::describe_posterior(bact_mods5b) # bio1 not sig
bayestestR::describe_posterior(bact_mods5c) # ai not sig 

# Compare the rsquared
brms::bayes_R2(bact_mods1) # 0.34
brms::bayes_R2(bact_mods2) # 0.27
brms::bayes_R2(bact_mods3) # 0.09
brms::bayes_R2(bact_mods4) # 0.33
brms::bayes_R2(bact_mods5a)# 0.22
brms::bayes_R2(bact_mods5b) # 0.22
brms::bayes_R2(bact_mods5c) # 0.22

# Getting deeper in algae models ----------------------------
# Compare models combining important variables both interacting and not interacting
# alg drought
# (fpom_mg_dry_log) (actual_water_mL_log) | (ai_log) | (bio1_log))   
alg_deep1 <- 
  brms::brm(formula = chl_a_ugL_log ~ 
              fpom_mg_dry_log + actual_water_mL_log + ai_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_drought, 
            family = gaussian())
alg_deep2 <- 
  update(alg_deep1,
         formula = chl_a_ugL_log ~
           fpom_mg_dry_log + actual_water_mL_log + bio1_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_drought, 
         family = gaussian())
alg_deep3 <- 
  update(alg_deep1,
         formula = chl_a_ugL_log ~
           (fpom_mg_dry_log + actual_water_mL_log)*ai_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_drought, 
         family = gaussian())
alg_deep4 <- 
  update(alg_deep1,
         formula = chl_a_ugL_log ~
           (fpom_mg_dry_log + actual_water_mL_log)*bio1_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_drought, 
         family = gaussian())

# alg fg
# (fpom_mg_dry_log) | actual_water_mL_log | bio1_log | bio15_log | ai_log | canopy_binary
## Take out biol15 because not present in most efficient parameters of drought expt data
## Keep canopy binary though
alg_deep5 <- 
  brms::brm(formula = alg_mL_log ~ 
              fpom_mg_dry_log + actual_water_mL_log + ai_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_fg, 
            family = gaussian())
alg_deep6 <- 
  update(alg_deep5,
         formula = alg_mL_log ~
           fpom_mg_dry_log + actual_water_mL_log + bio1_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
alg_deep7 <- 
  update(alg_deep5,
         formula = alg_mL_log ~
           (fpom_mg_dry_log + actual_water_mL_log)*ai_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())

# Takes a vrry long time to compute, forget
# alg_deep8 <- 
#   update(alg_deep5,
#          formula = alg_mL_log ~
#            (fpom_mg_dry_log + actual_water_mL_log)*bio1_log + (1 | site),
#          control = list(adapt_delta = 0.95,
#                         max_treedepth = 15),
#          iter = 5000,
#          newdata = equ_dats_fg, 
#          family = gaussian())

alg_deep9 <- 
  update(alg_deep5,
         formula = alg_mL_log ~
           fpom_mg_dry_log + actual_water_mL_log + canopy_binary + ai_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
alg_deep10 <- 
  update(alg_deep5,
         formula = alg_mL_log ~
           fpom_mg_dry_log + actual_water_mL_log + canopy_binary + bio1_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
alg_deep11 <- 
  update(alg_deep5,
         formula = alg_mL_log ~
           (fpom_mg_dry_log + actual_water_mL_log + canopy_binary)*ai_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())
alg_deep12 <- 
  update(alg_deep5,
         formula = alg_mL_log ~
           (fpom_mg_dry_log + actual_water_mL_log + canopy_binary)*bio1_log + (1 | site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_fg, 
         family = gaussian())

# LOO drought
## all equivalent
alg_deeploo1 <-
  list(alg_deep1, alg_deep2, alg_deep3) %>%
  purrr::map(~brms::loo(.x)) %>%
  brms::loo_compare() %>%
  get_eq_models()

# LOO FG
## all equivalent
alg_deeploo2 <-
  list(alg_deep5, alg_deep6, alg_deep7, alg_deep9,
       alg_deep10, alg_deep10, alg_deep11, alg_deep12) %>%
  purrr::map(~brms::loo(.x)) %>%
  brms::loo_compare() %>%
  get_eq_models()


# Getting deeper in bact models -------------------------------------------
# fpom_mg_dry_log +  alg_mL_log/chl_a_ugL_log | (ai_log) | (bio1_log) | (bio15_log)
# Let's keep the alg_mL because it has more data points
bact_deep1 <- 
  brms::brm(formula = bact_mL_log ~
              fpom_mg_dry_log + alg_mL_log + ai_log + (1 | site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = equ_dats_both, 
            family = gaussian())
bact_deep2 <- 
  update(bact_deep1,
         formula = bact_mL_log ~ 
           fpom_mg_dry_log + alg_mL_log + bio1_log + (1 | site),      
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_both, 
         family = gaussian())
bact_deep3 <- 
  update(bact_deep1,
         formula = bact_mL_log ~ 
           (fpom_mg_dry_log + alg_mL_log)*ai_log + (1 | site),      
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_both, 
         family = gaussian())
bact_deep4 <- 
  update(bact_deep1,
         formula = bact_mL_log ~ 
           (fpom_mg_dry_log + alg_mL_log)*bio1_log + (1 | site),      
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_both, 
         family = gaussian())
bact_deep5 <- 
  update(bact_deep1,
         formula = bact_mL_log ~ 
           fpom_mg_dry_log + bio1_log + (1 | site),      
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_both %>% 
           dplyr::filter(!is.na(alg_mL_log)), 
         family = gaussian())
# tried torun but very poor fit, let's forget about it
# bact_deep6 <- 
#   update(bact_deep1,
#          formula = bact_mL_log ~ 
#            fpom_mg_dry_log*bio1_log + (1 | site),      
#          control = list(adapt_delta = 0.95,
#                         max_treedepth = 15),
#          iter = 5000,
#          newdata = equ_dats_both %>% 
#            dplyr::filter(!is.na(alg_mL_log)), 
#          family = gaussian())
bact_deep6 <- 
  update(bact_deep1,
         formula = bact_mL_log ~ 
           fpom_mg_dry_log + ai_log + (1 | site),      
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_both %>% 
           dplyr::filter(!is.na(alg_mL_log)), 
         family = gaussian())
bact_deep7 <- 
  update(bact_deep1,
         formula = bact_mL_log ~ 
           fpom_mg_dry_log*ai_log + (1 | site),      
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = equ_dats_both %>% 
           dplyr::filter(!is.na(alg_mL_log)), 
         family = gaussian())

# LOO
## all equivalent
bact_deeploo <-
  list(bact_deep1, bact_deep2, bact_deep3, bact_deep4,
       bact_deep5, bact_deep6, bact_deep7) %>%
  purrr::map(~brms::loo(.x)) %>%
  brms::loo_compare() %>%
  get_eq_models()


# Use a priori judgement to pick best models ------------------------------
# Algae
## More data points and interactions between site and local variables
alg_final <- 
  brms::brm(formula = alg_mL_log ~
              (fpom_mg_dry_log + actual_water_mL_log)*ai_log + (1 | site),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 10000,
            data = equ_dats_fg, 
            family = gaussian())
brms::bayes_R2(alg_final) # 0.32

# Bacteria 
## Only one site variable, but still interaction
bact_final <- 
  brms::brm(formula = bact_mL_log ~
              fpom_mg_dry_log*ai_log + (1 | site),
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 10000,
            data = equ_dats_both, 
            family = gaussian())
brms::bayes_R2(bact_final) # 0.29

# Do predictions ---------------------------------------------------------
# Prepare bromeliad data
bromeliads <- 
  fw_database$bromeliads %>% 
  ## Keep important columns
  dplyr::select(visit_id, bromeliad_id, dataset_name, 
                total_detritus, actual_water) %>% 
  ## Add site coordinates and climate variables
  dplyr::left_join(fw_database$datasets %>% 
                     ## just keep lat long and site name
                     dplyr::select(site = name, lat, lon = lng) %>% 
                     ## Add env variables
                     add_env(raster_folder = here::here("chelsa")) %>% 
                     ## Just keep the ones we want
                     dplyr::select(dataset_name = site, ai),
                   by = "dataset_name")


# Use models to predict algae and bacteria
## Finish a few things in the data
bromeliads <- 
  bromeliads %>% 
  ## Rename
  dplyr::rename(site = dataset_name,
                fpom_mg_dry = total_detritus,
                actual_water_mL = actual_water) %>% 
  ## Log
  dplyr::mutate(dplyr::across(c(fpom_mg_dry, actual_water_mL,ai),
                              ~ifelse(.x == 0, 0, log(.x)),
                              .names = "{.col}_log"))

# Now predict algae and bacteria
## If "gaussian", sample new levels from the (multivariate) normal distribution 
## implied by the group-level standard deviations and correlations
bromeliads_predicted <- 
  bromeliads %>% 
  ## Algae
  dplyr::bind_cols(bromeliads %>%
                     ### Predictions
                     predict(object = alg_final,
                             allow_new_levels = TRUE, ### new random effects
                             sample_new_levels = "gaussian") %>%  
                     ### Convert to tibble and rename for quick stick
                     tibble::as_tibble() %>% 
                     dplyr::rename_with(~paste0("alg_mL_log_", stringr::str_to_lower(.x)))) %>% 
  ## Bacteria
  dplyr::bind_cols(bromeliads %>%
                     ### Predictions
                     predict(object = bact_final,
                             allow_new_levels = TRUE, ### new random effects
                             sample_new_levels = "gaussian") %>% 
                     ### Convert to tibble and rename for quick stick
                     tibble::as_tibble() %>% 
                     dplyr::rename_with(~paste0("bact_mL_log_", stringr::str_to_lower(.x))))

# Unlog predictors and convert to biomass using average mass of 
# algae (3.1e-11) and bacteria (Ducklow 2000, 1.6e-13g)
bromeliads_predicted <- 
  bromeliads_predicted %>% 
  dplyr::mutate(dplyr::across(dplyr::contains(c("bact")),
                              ~exp(.x) * 0.00000000000016)) %>% 
  dplyr::mutate(dplyr::across(dplyr::contains(c("alg")),
                              ~exp(.x) * 0.000000000031)) %>% 
  dplyr::rename_at(vars(dplyr::everything()),
                   ~stringr::str_replace(string = .,
                                         pattern = "_mL_log",
                                         replacement = "_g"))
# Save data
readr::write_csv(bromeliads_predicted,
                 here::here("data",
                            "alg_bact_predicted.csv"))
