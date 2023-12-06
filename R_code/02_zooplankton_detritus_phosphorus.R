# Get raw FPOM CPOM
library(tidyverse)
library(here)
library(brms)
library(bayestestR)


# Prepare protozoa data ---------------------------------------------------
# Read and clean env datar
proto_env <- 
  readr::read_csv(here::here("data",
                             "protozoa_enviro_data.csv")) %>% 
  ## Make a binary canopy column
  dplyr::mutate(canopy_binary = ifelse(canopy_cover < 0.5,
                                       0, 1)) %>% 
  ## Keep relevant columns
  dplyr::select(sampling, bromeliad, leaf, site, treatment, 
                canopy_binary, actual_water_mL = volume, total_detritus_g = litter)


# Read in the two data sets
proto <- 
  readr::read_csv(here::here("data",
                             "protozoa_data.csv")) %>% 
  ## Convert two names with two groups
  dplyr::rename("Flagellate_21" = Flagellate_21_Ciliate_13,
                "Flagellate_34" = Flagellate_34_Ciliate_23) %>% 
  ## Do sums at the scale of the different groups, rotifers being the only predators
  dplyr::mutate(abundance_preds = rowSums(dplyr::across(dplyr::contains(c("Rotifer")))),
                abundance_amoeba = rowSums(dplyr::across(dplyr::contains(c("Amoeba")))),
                abundance_diatom = rowSums(dplyr::across(dplyr::contains(c("Diatom")))),
                abundance_flagellate = rowSums(dplyr::across(dplyr::contains(c("Flagellate")))),
                abundance_ciliate = rowSums(dplyr::across(dplyr::contains(c("Ciliate|Vorticella"))))) %>% 
  ## Keep columns of interest
  dplyr::select(sampling, bromeliad, leaf, site, treatment, 
                dplyr::contains("abundance")) %>% 
  ## aaaand remove the two extra extra
  dplyr::select(-abundance_algae, - abundance1) %>% 
  ## join the environmental data
  dplyr::left_join(proto_env,
                   by = c("sampling", "bromeliad", "leaf",
                          "site", "treatment" )) %>% 
  ## Remove NA sites
  dplyr::mutate(site = ifelse(is.na(site),
                              "other", "site"))  %>% 
  ## Log
  dplyr::mutate(dplyr::across(c(actual_water_mL, total_detritus_g,
                                dplyr::contains("abundance")),
                              ~as.numeric(log(.+0.001)),
                              .names = "{.col}_log"))


# Do equations for each, with all combination of each variable
## Unfortunately only one site
# Make all possible combinations of variables
vars <- 
  ## Generate all possible combinations of variables
  tidyr::crossing(
    a = c(NA, "total_detritus_g_log"),
    b = c(NA, "actual_water_mL_log"),
    c = c(NA, "canopy_binary")) %>% 
  ## Combine with a + to get model specs
  tidyr::unite("a",
               a:c, 
               sep = " + ") %>% 
  ## Remove the NAS
  dplyr::mutate(a = stringr::str_remove_all(string = a,
                                            pattern = " [+] NA|NA [+] ")) %>% 
  dplyr::filter(!stringr::str_detect(string = a,
                                    pattern = "NA")) %>% 
  ## Add the random effect
  dplyr::mutate(a = paste(a, "+ (1|leaf/bromeliad/site)")) %>% 
  ## Pull to make vector
  dplyr::pull() %>% 
  ## Keep unique combinations
  unique()

    
    

# Models on mesozooplankton (predators) -----------------------------------
# Make list to store all models
meso_mods <- 
  list()
# change vars to include mesozooplankton as response variable
tempvars <- 
  paste0("abundance_preds_log ~ ", vars)

# Make loop that will make all models and save them in the list
for(i in 1:length(tempvars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(tempvars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = proto, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      meso_mods[[length(meso_mods) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(tempvars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = proto, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      meso_mods[[length(meso_mods) + 1]] <- 
        model
  }
  
}

# Get loo
meso_loo <- 
  meso_mods %>% 
  purrr::map(~loo(.x))

# Add names
names(meso_loo) <- 
  meso_mods %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()


## extract n of different points
npoints <- 
  meso_loo %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list based on number of data points
  temp <- 
    meso_loo %>% 
    purrr::keep(~any(nrow(.x$pointwise) == npoints[i]))
  ## If we have more than one model, keep the equivalent ones
  if(length(temp) > 1)
  {assign(
    paste0("meso_loo_", i),
    meso_loo %>% 
      purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
      brms::loo_compare() %>% 
      get_eq_models())} else {
        ## If we have only one model just keep it
        assign(
          paste0("meso_loo_", i),
          meso_loo %>% 
            purrr::keep(~any(nrow(.x$pointwise) == npoints[i])))
      }
}

# Get most parsimonious model in each category
meso_loo_1 # (actual_water_mL_log)
meso_loo_2 # (canopy_binary)


# Compute these models one by one and examine
## (fpom_mg_dry_log)
meso_mod1 <- 
  brms::brm(formula = abundance_preds_log ~ actual_water_mL_log + 
              (1 | leaf/bromeliad/site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = proto, 
            family = gaussian())
## (actual_water_mL_log)
meso_mod2 <- 
  update(meso_mod1,
         formula = abundance_preds_log ~ canopy_binary + 
           (1 | leaf/bromeliad/site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = proto, 
         family = gaussian())

# Compare the effect of independent variables
bayestestR::describe_posterior(meso_mod1) # actual water not sig
bayestestR::describe_posterior(meso_mod2) # canopy  sig

# Compare the rsquared
brms::bayes_R2(meso_mod1) # 0.07
brms::bayes_R2(meso_mod2) # 0.09



# Models on amoeba -----------------------------------
# Make list to store all models
amoe_mods <- 
  list()
# change vars to include amoezooplankton as response variable
tempvars <- 
  paste0("abundance_amoeba_log ~ ", vars)

# Make loop that will make all models and save them in the list
for(i in 1:length(tempvars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(tempvars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = proto, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      amoe_mods[[length(amoe_mods) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(tempvars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = proto, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      amoe_mods[[length(amoe_mods) + 1]] <- 
        model
  }
  
}

# Get loo
amoe_loo <- 
  amoe_mods %>% 
  purrr::map(~loo(.x))

# Add names
names(amoe_loo) <- 
  amoe_mods %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()


## extract n of different points
npoints <- 
  amoe_loo %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list based on number of data points
  temp <- 
    amoe_loo %>% 
    purrr::keep(~any(nrow(.x$pointwise) == npoints[i]))
  ## If we have more than one model, keep the equivalent ones
  if(length(temp) > 1)
  {assign(
    paste0("amoe_loo_", i),
    amoe_loo %>% 
      purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
      brms::loo_compare() %>% 
      get_eq_models())} else {
        ## If we have only one model just keep it
        assign(
          paste0("amoe_loo_", i),
          amoe_loo %>% 
            purrr::keep(~any(nrow(.x$pointwise) == npoints[i])))
      }
}


# Get most parsimonious model in each category
amoe_loo_1 # (actual_water_mL_log)
amoe_loo_2 # (canopy_binary)


# Compute these models one by one and examine
## (fpom_mg_dry_log)
amoe_mod1 <- 
  brms::brm(formula = abundance_amoeba_log ~ actual_water_mL_log + 
              (1 | leaf/bromeliad/site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = proto, 
            family = gaussian())
## (actual_water_mL_log)
amoe_mod2 <- 
  update(amoe_mod1,
         formula = abundance_amoeba_log ~ canopy_binary + 
           (1 | leaf/bromeliad/site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = proto, 
         family = gaussian())

# Compare the effect of independent variables
bayestestR::describe_posterior(amoe_mod1) # actual water not sig
bayestestR::describe_posterior(amoe_mod2) # canopy not sig

# Compare the rsquared
brms::bayes_R2(amoe_mod1) # 0.07
brms::bayes_R2(amoe_mod2) # 0.07s



# Models on diatoms -----------------------------------
# Make list to store all models
diat_mods <- 
  list()
# change vars to include diatzooplankton as response variable
tempvars <- 
  paste0("abundance_diatom_log ~ ", vars)

# Make loop that will make all models and save them in the list
for(i in 1:length(tempvars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(tempvars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = proto, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      diat_mods[[length(diat_mods) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(tempvars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = proto, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      diat_mods[[length(diat_mods) + 1]] <- 
        model
  }
  
}

# Get loo
diat_loo <- 
  diat_mods %>% 
  purrr::map(~loo(.x))

# Add names
names(diat_loo) <- 
  diat_mods %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()


## extract n of different points
npoints <- 
  diat_loo %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list based on number of data points
  temp <- 
    diat_loo %>% 
    purrr::keep(~any(nrow(.x$pointwise) == npoints[i]))
  ## If we have more than one model, keep the equivalent ones
  if(length(temp) > 1)
      {assign(
        paste0("diat_loo_", i),
        diat_loo %>% 
          purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
          brms::loo_compare() %>% 
          get_eq_models())} else {
            ## If we have only one model just keep it
            assign(
              paste0("diat_loo_", i),
              diat_loo %>% 
                purrr::keep(~any(nrow(.x$pointwise) == npoints[i])))
    }
}

# Get most parsimonious model in each category
diat_loo_1 # (actual_water_mL_log)
diat_loo_2 # (canopy_binary)


# Compute these models one by one and examine
## (fpom_mg_dry_log)
diat_mod1 <- 
  brms::brm(formula = abundance_diatom_log ~ actual_water_mL_log + 
              (1 | leaf/bromeliad/site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = proto, 
            family = gaussian())
## (actual_water_mL_log)
diat_mod2 <- 
  update(diat_mod1,
         formula = abundance_diatom_log ~ canopy_binary + 
           (1 | leaf/bromeliad/site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = proto, 
         family = gaussian())

# Compare the effect of independent variables
bayestestR::describe_posterior(diat_mod1) # actual water not sig
bayestestR::describe_posterior(diat_mod2) # canopy not sig

# Compare the rsquared
brms::bayes_R2(diat_mod1) # 0.11
brms::bayes_R2(diat_mod2) # 0.11




# Models on flagellates -----------------------------------
# Make list to store all models
flag_mods <- 
  list()
# change vars to include flagzooplankton as response variable
tempvars <- 
  paste0("abundance_flagellate_log ~ ", vars)

# Make loop that will make all models and save them in the list
for(i in 1:length(tempvars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(tempvars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = proto, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      flag_mods[[length(flag_mods) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(tempvars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = proto, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      flag_mods[[length(flag_mods) + 1]] <- 
        model
  }
  
}

# Get loo
flag_loo <- 
  flag_mods %>% 
  purrr::map(~loo(.x))

# Add names
names(flag_loo) <- 
  flag_mods %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()


## extract n of different points
npoints <- 
  flag_loo %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list based on number of data points
  temp <- 
    flag_loo %>% 
    purrr::keep(~any(nrow(.x$pointwise) == npoints[i]))
  ## If we have more than one model, keep the equivalent ones
  if(length(temp) > 1)
  {assign(
    paste0("flag_loo_", i),
    flag_loo %>% 
      purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
      brms::loo_compare() %>% 
      get_eq_models())} else {
        ## If we have only one model just keep it
        assign(
          paste0("flag_loo_", i),
          flag_loo %>% 
            purrr::keep(~any(nrow(.x$pointwise) == npoints[i])))
      }
}

# Get most parsimonious model in each category
flag_loo_1 # (actual_water_mL_log)
flag_loo_2 # (canopy_binary)


# Compute these models one by one and examine
## (fpom_mg_dry_log)
flag_mod1 <- 
  brms::brm(formula = abundance_flagellate_log ~ actual_water_mL_log + 
              (1 | leaf/bromeliad/site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = proto, 
            family = gaussian())
## (actual_water_mL_log)
flag_mod2 <- 
  update(flag_mod1,
         formula = abundance_flagellate_log ~ canopy_binary + 
           (1 | leaf/bromeliad/site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = proto, 
         family = gaussian())

# Compare the effect of independent variables
bayestestR::describe_posterior(flag_mod1) # actual water not sig
bayestestR::describe_posterior(flag_mod2) # canopy not sig

# Compare the rsquared
brms::bayes_R2(flag_mod1) # 0.02
brms::bayes_R2(flag_mod2) # 0.02




# Models on ciliates # very poor fit -----------------------------------
# Make list to store all models
cili_mods <- 
  list()
# change vars to include cilizooplankton as response variable
tempvars <- 
  paste0("abundance_ciliate_log ~ ", vars)

# Make loop that will make all models and save them in the list
for(i in 1:length(tempvars)){
  # Let's keep track of things
  print(i)
  print(vars[i])
  
  # If i = 1 compute model (i.e. compile in stan)
  if(i == 1){
    model1 <- 
      brms::brm(formula = brmsformula(tempvars[i]),      
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = proto, 
                family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model1) >1.05))
      cili_mods[[length(cili_mods) + 1]] <- 
        model1
    
  }
  # If i > 1 we just need to update
  # Compute model
  if(i > 1){
    model <- 
      update(model1,
             formula = brmsformula(tempvars[i]),
             control = list(adapt_delta = 0.95,
                            max_treedepth = 15),
             newdata = proto, 
             family = gaussian())
    
    # If model fit was good
    if(!any(brms::rhat(model) > 1.05))
      cili_mods[[length(cili_mods) + 1]] <- 
        model
  }
  
}

# Get loo
cili_loo <- 
  cili_mods %>% 
  purrr::map(~loo(.x))

# Add names
names(cili_loo) <- 
  cili_mods %>% 
  purrr::map(~.x$formula$formula[3]) %>% 
  purrr::list_flatten()


## extract n of different points
npoints <- 
  cili_loo %>% 
  purrr::map(~nrow(.x$pointwise)) %>% 
  purrr::list_c() %>% 
  unique()

## Now do loo compare on each set of model
for(i in 1:length(npoints)){
  # Filter model in list based on number of data points
  temp <- 
    cili_loo %>% 
    purrr::keep(~any(nrow(.x$pointwise) == npoints[i]))
  ## If we have more than one model, keep the equivalent ones
  if(length(temp) > 1)
  {assign(
    paste0("cili_loo_", i),
    cili_loo %>% 
      purrr::keep(~any(nrow(.x$pointwise) == npoints[i])) %>% 
      brms::loo_compare() %>% 
      get_eq_models())} else {
        ## If we have only one model just keep it
        assign(
          paste0("cili_loo_", i),
          cili_loo %>% 
            purrr::keep(~any(nrow(.x$pointwise) == npoints[i])))
      }
}

# Get most parsimonious model in each category
cili_loo_1 # (actual_water_mL_log)
cili_loo_2 # (canopy_binary)


# Compute these models one by one and examine
## (fpom_mg_dry_log)
cili_mod1 <- 
  brms::brm(formula = abundance_ciliate_log ~ actual_water_mL_log + 
              (1 | leaf/bromeliad/site),      
            control = list(adapt_delta = 0.95,
                           max_treedepth = 15),
            iter = 5000,
            data = proto, 
            family = gaussian())
## (actual_water_mL_log)
cili_mod2 <- 
  update(cili_mod1,
         formula = abundance_ciliate_log ~ canopy_binary + 
           (1 | leaf/bromeliad/site),
         control = list(adapt_delta = 0.95,
                        max_treedepth = 15),
         iter = 5000,
         newdata = proto, 
         family = gaussian())

# Compare the effect of independent variables
bayestestR::describe_posterior(cili_mod1) # actual water not sig
bayestestR::describe_posterior(cili_mod2) # canopy not sig

# Compare the rsquared
brms::bayes_R2(cili_mod1) # 0.02
brms::bayes_R2(cili_mod2) # 0.02


# Vinicius' data ----------------------------------------------------------
# Read data
vinizoo <-
  readr::read_csv(here::here("data",
                             "vini_zooplankton.csv")) %>%
  ## select columns
  dplyr::select(species = Species, V1:Al14) %>% 
  ## Remove those already in database
  dplyr::filter(stringr::str_detect(string = species,
                                    pattern = "Elpidium|Harpaticoida|Dero")) %>% 
  ## Make column if meso or micro
  dplyr::mutate(group = ifelse(stringr::str_detect(string = species,
                                                 pattern = "Tropocyclops|Ectocyclops"),
                              "micro", "meso")) %>% 
  ## Pivot to long format
  tidyr::pivot_longer(cols = where(is.numeric),
                      names_to = "bromeliad_id",
                      values_to = "abundance")

# Unique num

  
# Get average biomass of all groups, combine, and do models ---------------
fwdata$datasets %>% 
  dplyr::filter(num_bromeliads == length(unique(vinizoo$bromeliad_id)))
# Macae2008


