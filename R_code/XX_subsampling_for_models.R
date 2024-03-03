# Random sumbsamples of data to see the impact of food web size
library(tidyverse)
source(here::here("R_code",
                  "functions.R"))
library(brms)
library(DescTools)
library(here)
library(bayestestR)
library(ggplot2)
library(cowplot)

# Looking a bit at the data -----------------------------------------------
# First see how many bromeliads in each species pool
names(fluxes) %>% 
  tibble::as_tibble() %>% 
  tidyr::separate(value,
                  into = c("site", "n_broms")) %>% 
  dplyr::group_by(site) %>% 
  dplyr::distinct() %>% 
  dplyr::tally() %>% 
  dplyr::arrange(desc(n)) %>% 
  print(n = 30)

# Then see how many species in each species pool
site_matrices %>% 
  purrr::map(~row.names(.)) %>% 
  purrr::map(~.[!stringr::str_detect(string = .,
                                     pattern = "algae|bromeliad|FPOM|CPOM|zooplankton|microorganisms")]) %>% 
  purrr::map(~length(.)) %>% 
  tibble::as.tibble() %>% 
  tidyr::pivot_longer(dplyr::everything(),
                      names_to = "site",
                      values_to = "n") %>% 
  dplyr::arrange(desc(n)) %>% 
  print(n = 30)
  

# Then in each bromeliad
fluxes_long %>% 
  tidyr::pivot_longer(c(resource, consumer),
                      names_to = "dummy",
                      values_to = "species") %>% 
  dplyr::select(site, bromeliad_id, species) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(species %notin% basal) %>% 
  dplyr::group_by(site, bromeliad_id) %>% 
  dplyr::tally() %>% 
  dplyr::arrange(desc(n)) %>% 
  print(n = 30)

# Make mock data with pattern ----------------------------------------------------------
# 15 sites, 20 - 100 bromeliads, 15 - 40 species in each site, 5 - 20 species in each
# 15 sites
mock_sites <- 
  paste0("site", 
         seq(from = 1, 
             to = 15,
             by = 1))

# Sample 16 - 40 species in each site
mock_species <- 
  tibble::tibble(
    ## Sample number of species from uniform distribution
    site = rep(mock_sites, 
               round(runif(n= 15, min = 16, max = 40)))) %>% 
  ## Give each species a unique name based on row number
  dplyr::mutate(bromeliad_id = paste0("species_",
                                      seq(from = 1,
                                          to = nrow(.),
                                          by = 1))) %>% 
  ## Nest species per site (makes a column of vectors)
  tidyr::nest(.by = site, 
              .key = "specpool")

# Sample between 20 - 100 bromeliads for each of the 15 sites
mock_broms <- 
  tibble::tibble(
    ## Sample number of bromeliads from uniform distribution
    site = rep(mock_sites, 
               round(runif(n  = 15, min  = 20, max = 100)))) %>% 
  ## Give each bromeliad a unique name based on row number
  dplyr::mutate(bromeliad_id = paste0("brom_",
                                      seq(from = 1,
                                          to = nrow(.),
                                          by = 1))) %>% 
  ## Sample 5 - 15 species from the site pool for each bromeliad
  dplyr::rowwise() %>% 
  ## Do this by accessing the column of vectors for each site
  dplyr::mutate(community = list(sample(list_c(mock_species$specpool[which(mock_species$site == site)]) %>% 
                                         dplyr::pull(), 
                                       round(runif(n = 1, 
                                                   min = 5, 
                                                   max = 15)), 
                                       replace = FALSE, 
                                       prob = NULL))) %>% 
  ## Ungroup
  dplyr::ungroup()

# Number of bromeliads per site
mock_nbroms <- 
  mock_broms %>% 
  dplyr::select(site) %>% 
  dplyr::group_by(site) %>% 
  dplyr::tally()

# Number of species per site
mock_nspecs <- 
  mock_species %>% 
  tidyr::unnest_longer(specpool) %>% 
  dplyr::select(site) %>% 
  dplyr::group_by(site) %>% 
  dplyr::tally()

# Make data per consumer
## First unnest the data
mock_percons <- 
  mock_broms %>% 
  tidyr::unnest_longer(community) %>% 
  ## rename for ease of use
  dplyr::rename(consumer = community) %>% 
  ## Make columns that we will use later on, see distributions below 
  ### Each time add a random intercept based on site
  ### CPI crisp
  dplyr::mutate(cons_CPI_crisp = rnbinom(n = nrow(.),
                                         size = 1,
                                         mu = rep(round(runif(n = 15, 
                                                               min = 125, 
                                                               max = 175)), 
                                                  each = mock_nspecs$n[which(mock_nspecs$site == site)])) + 10) %>% 
  ### P content
  #### Scale from cons CPI
  dplyr::mutate(P_cons = (cons_CPI_crisp - min(cons_CPI_crisp))/(10*(max(cons_CPI_crisp) - min(cons_CPI_crisp))))  %>% 
  #### add noise
  dplyr::rowwise() %>% 
  dplyr::mutate(P_cons = P_cons + rnorm(n = 1,
                                        mean = 0,
                                        sd = 0.02)) %>% 
  dplyr::ungroup() %>% 
  #### just make sure it stays positive
  dplyr::mutate(P_cons = P_cons + abs(min(P_cons)) + 0.02) %>% 
  ### Raw fluxes
  #### Scale from P_cons
  dplyr::mutate(flux_raw = (P_cons - min(P_cons))/(max(P_cons) - min(P_cons)))  %>% 
  #### add noise
  dplyr::rowwise() %>% 
  dplyr::mutate(flux_raw = flux_raw + rnorm(n = 1, 
                                            mean = 0,
                                            sd = 0.1)) %>% 
  dplyr::ungroup() %>%  
  #### just make sure it stays positive
  dplyr::mutate(flux_raw = flux_raw + abs(min(flux_raw)) + 0.02) %>% 
  ## Fake C, N, P content
  dplyr::mutate(N_cons = runif(n = nrow(.), 
                          min = 0, 
                          max = 0.25),
                C_cons =  1 - P_cons - N_cons) %>% 
  ## Finally get flux proportin
  dplyr::group_by(site, bromeliad_id) %>% 
  dplyr::mutate(flux_prop = flux_raw/sum(flux_raw)) %>% 
  dplyr::ungroup()
 
  
  
  
# CPI crisp
hist(fluxes_long_percons$cons_CPI_crisp)
hist(rnbinom(n = 8000,
             size = 1,
             mu = 150))
hist(mock_percons$cons_CPI_crisp)
# P content
hist(fluxes_long_percons$P_cons)
hist(mock_percons$P_cons)
# Raw fluxes
hist(fluxes_long_percons$flux_raw)
hist(mock_percons$flux_prop)
# Relationships
plot(mock_percons$cons_CPI_crisp, mock_percons$P_cons)
plot(mock_percons$P_cons, mock_percons$flux_raw)



# Functions to return list of species to remove or add ----------------------------
# Remove
remove_specs <- function(data, number, mock_nspecs){
  
  
  # Select species in each
  remove <- 
    data %>% 
    dplyr::select(site, consumer) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(site) %>% 
    ## Either pick number given, or maximum number of species in species pool if lower
    slice_sample(n = ifelse(number < nrow(.), 
                            number, nrow(.)),
                 replace = F) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(consumer) %>% 
    dplyr::pull()
  
  # Remove species from data
  ret <- 
    data %>% 
    dplyr::filter(consumer %notin% remove)
  
  # Return
  return(ret)
}

# Add
add_specs <- function(data, number, mock_nspecs){
  
  
  # Select species in each
  keep <- 
    data %>% 
    dplyr::select(site, consumer) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(site) %>% 
    ## Either pick number given, or maximum number of species in species pool if lower
    slice_sample(n = ifelse(number < nrow(.), 
                            number, nrow(.)),
                 replace = F) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(consumer) %>% 

  
  # Remove species from data
  ret <- 
    data %>% 
    dplyr::filter(consumer %in% keep)
  
  # Return
  return(ret)
}


# Subset for site level - power analysis ---------------------------------------------------
# Full data
mock_perbrom_full <- 
  mock_percons %>% 
  dplyr::group_by(site, bromeliad_id) %>%
  ## Summarise
  dplyr::summarise(dplyr::across(flux_raw,
                                 ~ DescTools::Gini(.),
                                 .names = "gini"),
                   dplyr::across(flux_raw,
                                 ~ sum(.),
                                 .names = "totalflow_Jy")) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(
    mock_percons %>% 
      dplyr::select(site, bromeliad_id, C_cons , N_cons, P_cons) %>% 
      dplyr::distinct() %>% 
      dplyr::nest_by(site, bromeliad_id) %>%
      dplyr::summarise(detr_dist = average_dist(data)),
    by = c("site", "bromeliad_id"))

# Base models
## Gini - distance
mockfull_gini_detrdist <- 
  brms::brm(log(gini) ~
              log(detr_dist) + (1|site),
            iter = 5000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = mock_perbrom_full)

## Total flow - distance
mockfull_totalflow_detrdist <- 
  brms::brm(log(totalflow_Jy) ~
              poly(log(detr_dist), degree = 2) + (1|site),
            iter = 5000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = mock_perbrom_full)

## Emtpy tibble to store data
power_analysis <- 
  tibble::tibble()

# Really big looop
for(j in 1:100){
  
  # Remove zero to 36 species from each web and run models
  for(i in seq(from = 1, to = 36,  by = 3)){
    
    ## Just to keep track
    print(paste("rep", j, "minus", i, "species"))
    
    ## Make data
    dats <- 
      mock_percons %>% 
             ## Remove
             remove_specs(data = .,
                          number = i,
                          mock_nspecs = mock_nspecs) %>% 
             dplyr::group_by(site, bromeliad_id) %>%
             ## Summarise
             dplyr::summarise(dplyr::across(flux_raw,
                                            ~ DescTools::Gini(.),
                                            .names = "gini"),
                              dplyr::across(flux_raw,
                                            ~ sum(.),
                                            .names = "totalflow_Jy")) %>% 
             dplyr::ungroup() %>% 
             dplyr::left_join(
               mock_percons %>% 
                 dplyr::select(site, bromeliad_id, C_cons , N_cons, P_cons) %>% 
                 dplyr::distinct() %>% 
                 dplyr::nest_by(site, bromeliad_id) %>%
                 dplyr::summarise(detr_dist = average_dist(data)),
               by = c("site", "bromeliad_id"))
    
    ## Run gini detrdist model
    gini_detrdist_temp <- 
      update(object = mockfull_gini_detrdist,
             formula = log(gini) ~
               log(detr_dist) + (1|site),
             iter = 5000,
             family = gaussian(),      
             control = list(adapt_delta = 0.9,
                            max_treedepth = 10),
             newdata = dats)
    
    ## Run total flow detrdist model
    totalflow_detrdist_temp <- 
      update(object = mockfull_totalflow_detrdist,
             formula = log(totalflow_Jy) ~
               poly(log(detr_dist), degree = 2) + (1|site),
             iter = 5000,
             family = gaussian(),      
             control = list(adapt_delta = 0.9,
                            max_treedepth = 10),
             newdata = dats)
    
    ## Store results in dataframe
    power_analysis <- 
      power_analysis %>% 
      dplyr::bind_rows(tibble::tibble(rep = c(j, 
                                      j),
                              n_removed = c(i, 
                                            i),
                              model = c("gini_detrdist", 
                                        "totalflow_detrdist"),
                              rsquared = c(list(brms::bayes_R2(gini_detrdist_temp)), 
                                           list(brms::bayes_R2(totalflow_detrdist_temp))),
                              pd = c(list(bayestestR::pd(gini_detrdist_temp)),
                                     list(bayestestR::pd(totalflow_detrdist_temp))),
                              prcnt_rope = c(list(bayestestR::rope(gini_detrdist_temp)),
                                             list(bayestestR::rope(totalflow_detrdist_temp)))))
  
    
    
  }
  
  
  
    
    
}

# Save data
readr::write_rds(power_analysis,
                 here::here("data",
                            "power_analysis.rds"))



# Plot data from power analysis -------------------------------------------
# Read data
power_analysis_raw <- 
  readr::read_rds(here::here("data",
                            "power_analysis.rds"))

# Transform data to plot
power_analysis_plot <- 
  power_analysis_raw %>% 
    ## Values will change at each row
    dplyr::rowwise() %>% 
    ## Keep only the rsquared value
    dplyr::mutate(rsquared = rsquared[[1]]) %>% 
    ## Keep only the pd values
    dplyr::mutate(pd = list(pd$pd)) %>% 
    ## unnest the two or three values
    tidyr::unnest_wider(pd,
                        names_sep = c("_")) %>% 
    ## Rename the three columns
    dplyr::rename(pd_int = pd_1,
                  pd_rel1 = pd_2,
                  pd_rel2 = pd_3) %>% 
    ## Values will change at each row (disappeared somehow)
    dplyr::rowwise() %>% 
    ## Keep only the % in ROPE values
    dplyr::mutate(prcnt_rope = list(prcnt_rope$ROPE_Percentage)) %>% 
    ## unnest the two or three values
    tidyr::unnest_wider(prcnt_rope,
                        names_sep = c("_")) %>% 
    ## Rename the three columns
    dplyr::rename(rope_int = prcnt_rope_1,
                  rope_rel1 = prcnt_rope_2,
                  rope_rel2 = prcnt_rope_3)

# Make plots for gini vs. detritivore distance
## Rsquared
p1 <- 
  explo_line(dats = power_analysis_plot %>% 
                    dplyr::filter(model == "gini_detrdist"),
           x_var = "n_removed",
           y_var = "rsquared")
## pd of intercept
p2 <- 
  explo_line(dats = power_analysis_plot %>% 
               dplyr::filter(model == "gini_detrdist"),
             x_var = "n_removed",
             y_var = "pd_int")
## pd of linear term
p3 <- 
  explo_line(dats = power_analysis_plot %>% 
               dplyr::filter(model == "gini_detrdist"),
             x_var = "n_removed",
             y_var = "pd_rel1")
## rope of intercept
p4 <- 
  explo_line(dats = power_analysis_plot %>% 
               dplyr::filter(model == "gini_detrdist"),
             x_var = "n_removed",
             y_var = "rope_int")
## rope of linear term
p5 <- 
  explo_line(dats = power_analysis_plot %>% 
               dplyr::filter(model == "gini_detrdist"),
             x_var = "n_removed",
             y_var = "rope_rel1")
## Combine plot
### Little plot will be alone on top
row1 <- 
  cowplot::plot_grid(p1,
                     nrow = 1,
                     ncol = 1)
### Combine with rest








# Make mock data without pattern ------------------------------------------
# 15 sites, 100 - 500 bromeliads, 90 - 240 species in each site, 30 - 50 species in each
# 15 sites
mock_sites2 <- 
  paste0("site", 
         seq(from = 1, 
             to = 15,
             by = 1))

# Sample 90 - 240 species in each site
mock_species2 <- 
  tibble::tibble(
    ## Sample number of species from uniform distribution
    site = rep(mock_sites2, 
               round(runif(n= 15, min = 90, max = 240)))) %>% 
  ## Give each species a unique name based on row number
  dplyr::mutate(bromeliad_id = paste0("species_",
                                      seq(from = 1,
                                          to = nrow(.),
                                          by = 1))) %>% 
  ## Nest species per site (makes a column of vectors)
  tidyr::nest(.by = site, 
              .key = "specpool")

# Sample between 100 - 500 bromeliads for each of the 15 sites
mock_broms2 <- 
  tibble::tibble(
    ## Sample number of bromeliads from uniform distribution
    site = rep(mock_sites2, 
               round(runif(n  = 15, min  = 100, max = 500)))) %>% 
  ## Give each bromeliad a unique name based on row number
  dplyr::mutate(bromeliad_id = paste0("brom_",
                                      seq(from = 1,
                                          to = nrow(.),
                                          by = 1))) %>% 
  ## Sample 30 - 50 species from the site pool for each bromeliad
  dplyr::rowwise() %>% 
  ## Do this by accessing the column of vectors for each site
  dplyr::mutate(community = list(sample(list_c(mock_species2$specpool[which(mock_species2$site == site)]) %>% 
                                          dplyr::pull(), 
                                        round(runif(n = 1, 
                                                    min = 30, 
                                                    max = 50)), 
                                        replace = FALSE, 
                                        prob = NULL))) %>% 
  ## Ungroup
  dplyr::ungroup()

# Number of bromeliads per site
mock_nbroms2 <- 
  mock_broms2 %>% 
  dplyr::select(site) %>% 
  dplyr::group_by(site) %>% 
  dplyr::tally()

# Number of species per site
mock_nspecs2 <- 
  mock_species2 %>% 
  tidyr::unnest_longer(specpool) %>% 
  dplyr::select(site) %>% 
  dplyr::group_by(site) %>% 
  dplyr::tally()

# Make data per consumer
## First unnest the data
mock_percons2 <- 
  mock_broms2 %>% 
  tidyr::unnest_longer(community) %>% 
  ## rename for ease of use
  dplyr::rename(consumer = community) %>% 
  ## Make columns that we will use later on, see distributions below 
  ### Each time add a random intercept based on site
  ### CPI crisp
  dplyr::mutate(cons_CPI_crisp = rnbinom(n = nrow(.),
                                         size = 1,
                                         mu = 150)) %>% 
  ### P content
  dplyr::mutate(P_cons = runif(n = nrow(.),
                               min = 0.001,
                               max = 0.01)) %>% 
  ### Raw fluxes
  dplyr::mutate(flux_raw = rgamma(n = nrow(.),
                                  shape = 1,
                                  scale = 0.1))  %>% 
  ## Fake C, N, P content
  dplyr::mutate(N_cons = runif(n = nrow(.), 
                               min = 0, 
                               max = 0.25),
                C_cons =  1 - P_cons - N_cons) %>% 
  ## Finally get flux proportion
  dplyr::group_by(site, bromeliad_id) %>% 
  dplyr::mutate(flux_prop = flux_raw/sum(flux_raw)) %>% 
  dplyr::ungroup()




# CPI crisp
hist(fluxes_long_percons$cons_CPI_crisp)
hist(mock_percons2$cons_CPI_crisp)
# P content
hist(fluxes_long_percons$P_cons)
hist(mock_percons2$P_cons)
# Raw fluxes
hist(fluxes_long_percons$flux_raw)
hist(mock_percons2$flux_prop)
# Relationships
plot(mock_percons2$cons_CPI_crisp, mock_percons2$P_cons)
plot(mock_percons2$P_cons, mock_percons2$flux_raw)

# Increase at site level - artifacting ---------------------------------------------------
# Full data
mock_perbrom_full2 <- 
  mock_percons2 %>% 
  dplyr::group_by(site, bromeliad_id) %>%
  ## Summarise
  dplyr::summarise(dplyr::across(flux_raw,
                                 ~ DescTools::Gini(.),
                                 .names = "gini"),
                   dplyr::across(flux_raw,
                                 ~ sum(.),
                                 .names = "totalflow_Jy")) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(
    mock_percons2 %>% 
      dplyr::select(site, bromeliad_id, C_cons , N_cons, P_cons) %>% 
      dplyr::distinct() %>% 
      dplyr::nest_by(site, bromeliad_id) %>%
      dplyr::summarise(detr_dist = average_dist(data)),
    by = c("site", "bromeliad_id"))

# Base models
## Gini - distance
nopattern_gini_detrdist <- 
  brms::brm(log(gini) ~
              log(detr_dist) + (1|site),
            iter = 5000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = mock_perbrom_full2)

## Total flow - distance
nopattern_totalflow_detrdist <- 
  brms::brm(log(totalflow_Jy) ~
              poly(log(detr_dist), degree = 2) + (1|site),
            iter = 5000,
            family = gaussian(),      
            control = list(adapt_delta = 0.9,
                           max_treedepth = 10),
            data = mock_perbrom_full2)

## Emtpy tibble to store data
artifacting <- 
  tibble::tibble()

# Really big looop
for(j in 1:100){
  
  # Remove zero to 36 species from each web and run models
  for(i in seq(from = 1, to = 70,  by = 3)){
    
    ## Just to keep track
    print(paste("rep", j, "with", i, "species"))
    
    ## Make data
    dats <- 
      mock_percons2 %>% 
      ## Remove
      add_specs(data = .,
                   number = i,
                   mock_nspecs = mock_nspecs2) %>% 
      dplyr::group_by(site, bromeliad_id) %>%
      ## Summarise
      dplyr::summarise(dplyr::across(flux_raw,
                                     ~ DescTools::Gini(.),
                                     .names = "gini"),
                       dplyr::across(flux_raw,
                                     ~ sum(.),
                                     .names = "totalflow_Jy")) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(
        mock_percons2 %>% 
          dplyr::select(site, bromeliad_id, C_cons , N_cons, P_cons) %>% 
          dplyr::distinct() %>% 
          dplyr::nest_by(site, bromeliad_id) %>%
          dplyr::summarise(detr_dist = average_dist(data)),
        by = c("site", "bromeliad_id"))
    
    ## Run gini detrdist model
    gini_detrdist_temp <- 
      update(object = mockfull_gini_detrdist,
             formula = log(gini) ~
               log(detr_dist) + (1|site),
             iter = 5000,
             family = gaussian(),      
             control = list(adapt_delta = 0.9,
                            max_treedepth = 10),
             newdata = dats)
    
    ## Run total flow detrdist model
    totalflow_detrdist_temp <- 
      update(object = mockfull_totalflow_detrdist,
             formula = log(totalflow_Jy) ~
               poly(log(detr_dist), degree = 2) + (1|site),
             iter = 5000,
             family = gaussian(),      
             control = list(adapt_delta = 0.9,
                            max_treedepth = 10),
             newdata = dats)
    
    ## Store results in dataframe
    artifacting <- 
      artifacting %>% 
      dplyr::bind_rows(tibble::tibble(rep = c(j, 
                                              j),
                                      n_removed = c(i, 
                                                    i),
                                      model = c("gini_detrdist", 
                                                "totalflow_detrdist"),
                                      rsquared = c(list(brms::bayes_R2(gini_detrdist_temp)), 
                                                   list(brms::bayes_R2(totalflow_detrdist_temp))),
                                      pd = c(list(bayestestR::pd(gini_detrdist_temp)),
                                             list(bayestestR::pd(totalflow_detrdist_temp))),
                                      prcnt_rope = c(list(bayestestR::rope(gini_detrdist_temp)),
                                                     list(bayestestR::rope(totalflow_detrdist_temp)))))
    
    
    
  }
}

# Save data
readr::write_rds(artifacting,
                 here::here("data",
                            "artifacting.rds"))








