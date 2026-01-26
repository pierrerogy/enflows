# Make plots

# Load libraries 
library(tidyverse)
library(here)
library(cowplot)
source(here::here("R_code",
                  "functions.R"))


# Load data and models (ran on the cluster) -------------------------------
# Load data
fluxes_long<- 
  readr::read_csv(here::here("data",
                             "fluxes_long.csv"))
fluxes_long_percons <- 
  readr::read_csv(here::here("data",
                             "fluxes_long_percons.csv"))
fluxes_long_perpair <-
  readr::read_csv(here::here("data",
                             "fluxes_long_perpair.csv"))  
fluxes_long_perbrom <-
  readr::read_csv(here::here("data",
                             "fluxes_long_perbrom.csv"))

# Load models
consflow_model <- 
  readr::read_rds(here::here("data",
                             "consflow_model.RDS"))
cpi_model <- 
  readr::read_rds(here::here("data",
                             "cpi_model.RDS"))
pairflow_model <- 
  readr::read_rds(here::here("data",
                             "pairflow_model.RDS"))
totalflow_model_vol <- 
  readr::read_rds(here::here("data",
                             "totalflow_model_vol.RDS"))
totalflow_model_dist <- 
  readr::read_rds(here::here("data",
                             "totalflow_model_dist.RDS"))
predflow_model_vol <- 
  readr::read_rds(here::here("data",
                             "predflow_model_vol.RDS"))
predflow_model_dist <- 
  readr::read_rds(here::here("data",
                             "predflow_model_dist.RDS"))
gini_model_vol <- 
  readr::read_rds(here::here("data",
                             "gini_model_vol.RDS"))
gini_model_dist <- 
  readr::read_rds(here::here("data",
                             "gini_model_dist.RDS"))





# Random plotting of some food webs -----------------------------------------------------
# Plot whatever we want
plot_food_web(bromtoplot = "7606",
              sitetoplot = "PICIN",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)


# Plot most unequal bromeliad
## Find which ones are in the top 10%
fluxes_long_perbrom %>% 
  dplyr::slice_max(order_by = gini,
                   prop = 0.1) %>% 
  print(n = 150)
## Plot
plot_food_web(bromtoplot = "8021",
              sitetoplot = "GAMAS",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)

# Plot least unequal bromeliads
## Find which ones are in the bottom 10%
fluxes_long_perbrom %>% 
  dplyr::slice_min(order_by = gini,
                   prop = 0.1) %>% 
  print(n = 150)
## Plot
plot_food_web(bromtoplot = "5181",
              sitetoplot = "PTL",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)

# Plot something in the middle
## Find closest one to the mean
fluxes_long_perbrom %>% 
  dplyr::filter(gini == DescTools::Closest(gini, 
                                           mean(fluxes_long_perbrom$gini,
                                                na.rm = TRUE),
                                           na.rm = TRUE))
## Plot
plot_food_web(bromtoplot = "7306",
              sitetoplot = "SISGA",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)


# Find some nice one in Pitilla
## List of ids
fluxes_long %>% 
  dplyr::filter(site_code == "PTL") %>% 
  dplyr::select(bromeliad_id) %>% 
  dplyr::distinct() %>%  
  dplyr::pull()
## Plot
plot_food_web(bromtoplot = "346",
              sitetoplot = "PTL",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)
## let's select 346

# Plotting PTL 346-----------------------------------------------------
# Get the data
web <- 
  fluxes_long %>% 
  dplyr::filter(site_code == "PTL" & bromeliad_id == "346") %>% 
  ## Remove cannibalism
  dplyr::filter(consumer != resource)

# Define basal species
basal <- 
  c("algae", "FPOM", "CPOM", "microzooplankton",
    "mesozooplankton", "bromeliad.tissue", "microorganisms")

# Get list of species in the web
speclist <- 
  unique(c(web$consumer,
           web$resource))

# Get site interaction matrix
site_mat <- 
  binary_matrices[[which(names(binary_matrices) == "PTL")]]

# Get trophic levels within this bromeliad
levels <- 
  ## Subset site matrix by what is inside the bromeliad
  site_mat[rownames(site_mat) %in% speclist,
           colnames(site_mat) %in% speclist] %>% 
  ## Use function to get trophic levels
  TL() %>% 
  tibble::as_tibble(rownames = "species") %>% 
  dplyr::rename(TL = V1) %>%
  # Round trophic level to nearest 0.5
  dplyr::mutate(TL  = round((TL * 2)) / 2) %>% 
  # Modify that for micro- and mesozooplankton %>% 
  dplyr::mutate(TL = ifelse(stringr::str_detect(string = species,
                                                pattern = "plankton"),
                            1.5, TL)) %>% 
  ## Give random x axis position for organisms
  ## Assign x value to each species based on trophic level
  dplyr::mutate(pos = get_TL_pos(dats = .)) %>%
   ## Add biomass
  dplyr::left_join(clean_biomass_data %>% 
                     dplyr::filter(site_code == "PTL" & bromeliad_id == "346") %>% 
                     dplyr::select(biomass = tot_biomass_g,
                                   species = bwg_name),
                   by = "species") %>% 
  ## Add label, done inefficiently
  dplyr::mutate(label = c("Algae", "FPOM", "CPOM", "Bacteria",
                          "Micro", "Meso","Scirtidae", "italic(Polypedilum)",
                          "Orthocladiinae", "italic(W.~melanopus)", 
                          "italic(C.~jenningsi)", "italic(Copestylum)",
                          "italic(A.~kompii)", "italic(A.~neivai)", 
                          "italic(Trentepohlia)", "Periscelididae",
                          "italic(Monopelopia)", "italic(S.~fulvohirtum)", 
                          "italic(M.~modesta)"))

# Add trophic level and biomass to web
web <- 
  web %>% 
  ## For resources
  dplyr::left_join(levels %>% 
                     dplyr::rename(res_TL = TL,
                                   resource = species,
                                   res_pos = pos,
                                   res_biomass = biomass),
                   by = "resource") %>% 
  ## For consumers
  dplyr::left_join(levels %>% 
                     dplyr::rename(cons_TL = TL,
                                   consumer = species,
                                   cons_pos = pos,
                                   cons_biomass = biomass),
                   by = "consumer")

# Plot
## Base plot on levels tibble
pitiweb <- 
  ggplot(data = levels, 
         aes(y  = TL, 
             x = pos)) +
  ## Add curves between consumer resource pairs
  ## with shading and colour based on strength of flux
  geom_curve(data = web %>% 
               dplyr::filter(!is.na(cons_biomass)), 
             aes(y = res_TL, yend = cons_TL,
                 x = res_pos, xend = cons_pos,
                 colour = flux_raw, alpha = flux_raw),
             curvature = 0.1, arrow = arrow(length = unit(0.01, "npc"),
                                            ends = "last", type = "closed"),
             lineend = 'butt', linewidth = 1.5) +
  ## Give better colours
  scale_colour_viridis(option = "inferno",
                       end = 0.8,
                       name = expression(paste("Energy \nflux (J.y"^"-1"*")"))) +
  ## Remove shading legend
  scale_alpha_continuous(guide = 'none') +
  ## Axis labels
  xlab("") +
  ylab("Trophic level") +
  ## Species labels
  geom_label(aes(label = label),
             color = "black",
             fill = "white",
             size = 4,
             parse = TRUE) +
  xlim(-1.1, 1.1)+
  ## Styling
  theme(legend.justification = NULL,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))

## Save 
ggplot2::ggsave(here::here("plots",
                           "fig1.pdf"),
                pitiweb,
                height = 10,
                width = 10)
ggplot2::ggsave(here::here("plots",
                           "fig1.jpeg"),
                pitiweb,
                height = 10,
                width = 10)

# Figure 2 (run script 04) ------------------------------------------------
# Fig 2a
## Plot
fig2a <- 
  basic_line(model = consflow_model, 
             index = 1,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-3), 
             x_var = "cons_CPI_defuzzified", 
             y_var= "flux_raw") +
  scale_x_continuous(breaks = c(3, 4, 5, 6),
                     labels= c(20, 50, 150, 400)) +
  scale_y_continuous(limits = c(-10, 10),
                     breaks = seq(from = -10, 
                                  to = 10,
                                  by = 5),
                     labels= c(expression(paste("10"^{-5})),
                               0.001,
                               1,
                               150,
                               2000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))

# Fig 2b
## Plot
fig2b <- 
  basic_line(model = cpi_model, 
             index = 1,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-3), 
             x_var = "avg_bodymass_g", 
             y_var= "cons_CPI_defuzzified") +
  scale_x_continuous(limits = c(-15, 0),
                     breaks = seq(from = -15, 
                                  to = 0,
                                  by = 5),
                     labels= c(expression(paste("10"^{-5})),
                               expression(paste("10"^{-5})),
                               0.006,
                               1)) +
  scale_y_continuous(breaks = c(3, 4, 5, 6),
                     labels= c(20, 50, 150, 400)) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))
 
# Compile
fig2 <- 
  cowplot::plot_grid(plotlist = list(fig2a +
                                       ggtitle("(a)"), 
                                     fig2b +
                                       ggtitle("(b)")),
                     nrow = 1)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig2.pdf"),
                fig2,
                height = 6,
                width = 10)

  



# Figure 3 (run script 04) ------------------------------------------------
# Plot
fig3 <- 
  basic_line(model = pairflow_model, 
             index = 1,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "func_dist",
             y_var= "flux_raw_pair",
             scale = "x_original") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 5)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 7)),
        axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.4)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig3.pdf"),
                fig3,
                height = 10,
                width = 10)

# Figure 4 (run script 04)  ------------------------------------------------
# Fig 4a
fig4a <- 
  basic_line(model = totalflow_model_vol, 
              index = 1,
              dats = fluxes_long_perbrom, 
              x_var = "prop_hull",
              y_var= "totalflow_Jy") +
  scale_y_continuous(breaks = c(-2.5,
                                0,
                                2.5,
                                5,
                                7),
                     labels = c(0.1,
                                1,
                                10,
                                100,
                                1000),
                     limits = c(-2.5, 7.2)) +
  scale_x_continuous(breaks = c(seq(from = -15,
                                    to = 0,
                                    by = 3)),
                     labels = c(expression(paste("10"^{-7})),
                                expression(paste("10"^{-6})),
                                expression(paste("10"^{-4})),
                                0.002,
                                0.05,
                                1)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4b
fig4b <- 
  basic_line(model = predflow_model_vol, 
             index = 1,
             dats = fluxes_long_perbrom %>% 
               dplyr::mutate(predflow_Jy = predflow_Jy + 0.001), 
             x_var = "prop_hull",
             y_var= "predflow_Jy") +
  scale_y_continuous(breaks = c(-6.9,
                                -5,
                                -2.5,
                                0,
                                2.5, 
                                5,
                                7),
                     labels = c(0,
                                0.005,
                                0.1,
                                1,
                                10,
                                150,
                                1000),
                     limits = c(-9.8, 7.2)) +
  scale_x_continuous(breaks = c(seq(from = -15,
                                    to = 0,
                                    by = 3)),
                     labels = c(expression(paste("10"^{-7})),
                                expression(paste("10"^{-6})),
                                expression(paste("10"^{-4})),
                                0.002,
                                0.05,
                                1)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4c
fig4c <- 
  basic_line(model = gini_model_vol, 
             index = 1,
             dats = fluxes_long_perbrom, 
             x_var = "prop_hull",
             y_var = "gini",
             scale = "y_original") +
  scale_x_continuous(breaks = c(seq(from = -15,
                                    to = 0,
                                    by = 3)),
                     labels = c(expression(paste("10"^{-7})),
                                expression(paste("10"^{-6})),
                                expression(paste("10"^{-4})),
                                0.002,
                                0.05,
                                1)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4d
fig4d <- 
  basic_line(model = totalflow_model_dist, 
             index = 1,
             dats = fluxes_long_perbrom, 
             x_var = "avg_func_dist",
             y_var= "totalflow_Jy",
             scale = "x_original") +
  scale_y_continuous(breaks = c(-2.5,
                                0,
                                2.5,
                                5,
                                7),
                     labels = c(0.1,
                                1,
                                10,
                                100,
                                1000),
                     limits = c(-2.5, 7.2)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4e
fig4e <- 
  basic_line(model = predflow_model_dist, 
             index = 1,
             dats = fluxes_long_perbrom %>% 
               dplyr::mutate(predflow_Jy = predflow_Jy + 0.001), 
             x_var = "avg_func_dist",
             y_var= "predflow_Jy",
             scale = "x_original") +
  scale_y_continuous(breaks = c(-6.9,
                                -5,
                                -2.5,
                                0,
                                2.5, 
                                5,
                                7),
                     labels = c(0,
                                0.005,
                                0.1,
                                1,
                                10,
                                150,
                                1000),
                     limits = c(-7.2, 7.2)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4f
fig4f <- 
  basic_line(model = gini_model_dist, 
             index = 1,
             dats = fluxes_long_perbrom, 
             x_var = "avg_func_dist",
             y_var= "gini",
             scale = "both_original") +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 


# Compile
fig4 <- 
  cowplot::plot_grid(plotlist = list(fig4a +
                                       ggtitle("(a)"), 
                                     fig4b +
                                       ggtitle("(b)"), 
                                     fig4c +
                                       ggtitle("(c)"), 
                                     fig4d +
                                       ggtitle("(d)"), 
                                     fig4e +
                                       ggtitle("(e)"), 
                                     fig4f +
                                       ggtitle("(f)")),
                     nrow = 2)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig4.pdf"),
                fig4,
                height = 10,
                width = 15)



# Figures S1 - s10 (run script 04)----------------------------------------------
# Load site matrices
site_matrices <- 
  load_matrices(folder = here::here("data",
                                    "matrices"),
                ending = "_PM.txt") 

# Load models
masterlist_sites <- 
  readr::read_rds(here::here("data",
                             "masterlist_sites.RDS"))

# Lists to store plots
masterlist_siteplots <- 
  list()

# Make list for x and y labels
## Xlabs
xlabs <- 
  c("cons_P",
    "CN_dist",
    "CP_dist",
    "CNP_dist",
    "avg_CN_dist",
    "avg_CP_dist",
    "avg_CNP_dist",
    "avg_CN_dist",
    "avg_CP_dist",
    "avg_CNP_dist")
## Ylabs
ylabs <- 
  c("flux_raw",
    "flux_raw_pair",
    "flux_raw_pair",
    "flux_raw_pair",
    "totalflow_Jy",
    "totalflow_Jy",
    "totalflow_Jy",
    "predflow_Jy",
    "predflow_Jy",
    "predflow_Jy")


# Loop to plot for all sites

for(i in 1:length(masterlist_sites)){
  ## Extract list
  sublist <- 
    masterlist_sites[[i]]
  
  ## Temporary list to store plots
  templist <- 
    list()
  
  ## Sanity check
  print(i)
  
  ## Prep data for plots
  if(i == 1){
    dats <- 
      fluxes_long_percons %>% 
      dplyr::filter(site == names(site_matrices)[i]) %>% 
      dplyr::mutate(flux_prop = flux_prop + 1e-3,
                    flux_raw = flux_raw + 1e-3)
  } else if(i %in% c(2:4)){
    dats <- 
      fluxes_long_perpair %>%
      dplyr::filter(site == names(site_matrices)[i]) %>%
      dplyr::mutate(flux_prop_pair = flux_prop_pair + 1e-3,
                    flux_raw_pair = flux_raw_pair + 1e-3)
  } else {
    dats <- 
      dats <- 
      fluxes_long_perbrom %>% 
      dplyr::filter(site == names(site_matrices)[i])
  }
  
  
  ## Loop running thRrough list of models
  for(j in 1:length(site_matrices)){
    ### Sanity check
    print(names(site_matrices)[j])
    
    ### Extract model
    model <- 
      sublist[[j]]
    
    ## Make title 
    title <- 
      paste0(names(site_matrices)[j],
             "; N = ",
             nrow(model$data))
    
    ## Make different rules based on index
    ### Here we have two relationships in the first model
    if(i == 1){
      tempplot <- 
        list(basic_line(model = model, 
                   dats = model$data, 
                   index = 1,
                   x_var = "cons_CPI_crisp", 
                   y_var = ylabs[i]) +
               ggtitle(title) +
               annotate("text",
                        x = -Inf, 
                        y = Inf,
                        hjust = 0,
                        vjust = 2,
                        label = title_sig(sitename = names(site_matrices)[j], 
                                          model = model, 
                                          index = 1)) +
               theme(axis.title.y = element_text(margin = margin(r = 0,  l = 10))),
             basic_line(model = model, 
                        index = 2,
                        dats = model$data, 
                        x_var = "cons_P", 
                        y_var = ylabs[i]) +
               ggtitle(title) +
               annotate("text",
                        x = -Inf, 
                        y = Inf,
                        hjust = 0,
                        vjust = 2,
                        label = title_sig(sitename = names(site_matrices)[j], 
                                          model = model, 
                                          index = 2)) +
               theme(axis.title.y = element_text(margin = margin(r = 0,  l = 10))))
    } else 
      ### Here we have the interaction
      if(i %in% c(2:3)){
        tempplot <- 
          basic_line(model = model, 
                     dats = model$data, 
                     index = 3,
                     colour = paste0("lower_cons_",
                                     stringr::str_replace(string = xlabs[i],
                                                   pattern = "_dist",
                                                   replacement = "")),
                     x_var = xlabs[i],
                     y_var = ylabs[i],
                     scale = "none") +
          ggtitle(title) + 
          theme(axis.title.y = element_text(margin = margin(r = 0,  l = 10),
                                            size = rel(0.9)),
                legend.position = "none")
        
      } else{
        ### Here we are all good
        #### Make plot
        tempplot <- 
          basic_line(model = model, 
                     dats = model$data,
                     index = 1,
                     x_var = xlabs[i], 
                     y_var = ylabs[i],
                     scale = "none") +
          ggtitle(title) + 
          annotate("text",
                   x = -Inf, 
                   y = Inf,
                   hjust = 0,
                   vjust = 2, 
                   label = title_sig(sitename = names(site_matrices)[j], 
                                     model = model, 
                                     index = 1)) +
          theme(axis.title.y = element_text(margin = margin(r = 0,  l = 15),
                                            size = rel(0.8)))

      }
  
    ## Add to list
    templist[[length(templist) + 1]] <- 
      tempplot
    
  }
  
  ## Add to list
  masterlist_siteplots[[length(masterlist_siteplots) + 1]] <- 
    templist
  
}

# Save plots
## fluxraw vs cpi
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[1]] %>% 
                       purrr::map(1),
                     nrow = 6,
                     ncol = 5) 
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs1.pdf")) 
## fluxraw vs p
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[1]] %>% 
                       purrr::map(2),
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs2.pdf"))
## fluxraw vs cndist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[2]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs3.pdf"))
## fluxraw vs cpdist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[3]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs4.pdf"))
## fluxraw vs cnpdist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[4]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs5.pdf"))
## totalflow vs cn dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[5]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs6.pdf"))
## totalflow vs cp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[6]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs7.pdf"))
## totalflow vs cnp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[7]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs8.pdf"))
## pred flow vs cn dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[8]],
                     nrow = 6,
                     ncol = 5,
                     align = "hv")
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs9.pdf"))
## pred flow vs cp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[9]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs10.pdf"))
## pred flow vs cnp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[10]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs11.pdf"))





