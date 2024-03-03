# Nice plot of a food web

# Load libraries 
library(tidyverse)
library(here)
library(cowplot)
source(here::here("R_code",
                  "functions.R"))

# Random plotting of some food webs -----------------------------------------------------
# Plot most unequal bromeliad
## Find which one it is
fluxes_long_perbrom %>% 
  dplyr::filter(gini == max(gini))
## Plot
plot_food_web(bromtoplot = "9151",
              sitetoplot = "SONLOW",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)

# Plot least unequal bromeliad
## Find which one it is
fluxes_long_perbrom %>% 
  dplyr::filter(gini == min(gini))
## Plot
plot_food_web(bromtoplot = "8001",
              sitetoplot = "GAMAS",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)

# Plot something in the middle
## Find closest one to the mean
fluxes_long_perbrom %>% 
  dplyr::filter(gini == DescTools::Closest(gini, 
                                           mean(fluxes_long_perbrom$gini)))
## Plot
plot_food_web(bromtoplot = "12006",
              sitetoplot = "MAC",
              fluxes_long_format = fluxes_long,
              clean_biomass_data = biomass_clean,
              binary_matrices = site_matrices)


# Find some nice one in Pitilla
## List of ids
fluxes_long %>% 
  dplyr::filter(site == "PTL") %>% 
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
  dplyr::filter(site == "PTL" & bromeliad_id == "346") %>% 
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









# Figure 2 (run script 03) ------------------------------------------------
## Plot by hand because here covariate
# Get effect
model_effect <- 
  brms::conditional_effects(cpi_model_P,
                            method = "fitted")$cons_P 
fig2a <- 
  ggplot(data = fluxes_long_percons,
         aes(x = jitter(log(cons_P)),
             y = log(cons_CPI_crisp))) +
  ## Add raw data points
  geom_jitter(colour = "lightsteelblue4",
              alpha = 0.4) +
  ## Line
  geom_line(data = model_effect,
            aes(x = log(effect1__),
                y = estimate__),
            linewidth = 1.2,
            colour = "steelblue") +
  ## Confidence interval
  geom_ribbon(data = model_effect,
              aes(x = log(cons_P),
                  ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ## Styling
  xlab(full_label("cons_P")) +
  ylab(full_label("cons_CPI_crisp")) + 
  scale_y_continuous(breaks = c(3, 4, 5, 6),
                     labels= c(20, 50, 150, 400)) +
  scale_x_continuous(breaks = seq(from = -5.5, 
                                  to = -3.5,
                                  by = 0.5),
                     labels= round(exp(seq(from = -5.5, 
                                           to = -3.5,
                                           by = 0.5)),
                                   digits = 3)) +

  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

# Fig 2b
fig2b <- 
  basic_line(model = fluxprop_model_cpi, 
             dats = fluxes_long_percons, 
             x_var = "cons_CPI_crisp", 
             y_var= "flux_prop") +
  scale_y_continuous(limits = c(-20, 1),
                     breaks = seq(from = -20, 
                                  to = 0,
                                  by = 5),
                     labels= c(0,
                               expression(paste("10"^{-7})),
                               expression(paste("10"^{-5})),
                               0.001,
                               1)) + 
  scale_x_continuous(breaks = c(3, 4, 5, 6),
                     labels= c(20, 50, 150, 400))


# Fig 2c
fig2c <- 
  basic_line(model = fluxprop_model_p, 
             dats = fluxes_long_percons, 
             x_var = "cons_P", 
             y_var = "flux_prop") +
  scale_y_continuous(limits = c(-20, 1),
                     breaks = seq(from = -20, 
                                  to = 0,
                                  by = 5),
                     labels= c(0,
                               expression(paste("10"^{-7})),
                               expression(paste("10"^{-5})),
                               0.001,
                               1)) +
  scale_x_continuous(breaks = seq(from = -5.5, 
                                  to = -3.5,
                                  by = 0.5),
                     labels= round(exp(seq(from = -5.5, 
                                           to = -3.5,
                                           by = 0.5)),
                                   digits = 3))

# Fig 2d
fig2d <- 
  basic_line(model = fluxraw_model_cpi, 
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 0.0000000001), 
             x_var = "cons_CPI_crisp", 
             y_var= "flux_raw") + 
  scale_x_continuous(breaks = c(3, 4, 5, 6),
                     labels= c(20, 50, 150, 400)) +
  scale_y_continuous(limits = c(-15, 10),
                     breaks = seq(from = -15, 
                                  to = 10,
                                  by = 5),
                     labels= c(expression(paste("10"^{-7})),
                               expression(paste("10"^{-5})),
                               0.001,
                               1,
                               150,
                               2000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))

# Fig 2e
fig2e <- 
  basic_line(model = fluxraw_model_p, 
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 0.0000000001), 
             x_var = "cons_P", 
             y_var = "flux_raw",
             scale = "log") +
  scale_y_continuous(limits = c(-15, 10),
                     breaks = seq(from = -15, 
                                  to = 10,
                                  by = 5),
                     labels= c(expression(paste("10"^{-7})),
                               expression(paste("10"^{-5})),
                               0.001,
                               1,
                               150,
                               2000)) +
  scale_x_continuous(breaks = seq(from = -5.5, 
                                  to = -3.5,
                                  by = 0.5),
                     labels= round(exp(seq(from = -5.5, 
                                           to = -3.5,
                                           by = 0.5)),
                                   digits = 3)) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))
 
# Compile
fig2 <- 
  cowplot::plot_grid(plotlist = list(fig2a +
                                       ggtitle("(a)"), 
                                     fig2b +
                                       ggtitle("(b)"), 
                                     fig2c +
                                       ggtitle("(c)"),
                                     fig2d +
                                       ggtitle("(d)"),
                                     fig2e +
                                       ggtitle("(e)")),
                     nrow = 3)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig2.pdf"),
                fig2,
                height = 15,
                width = 10)

  



# Figure 3 (run script 03) ------------------------------------------------
# Fig 3a
fig3a <- 
  basic_line(model = fluxprop_model_cndist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_prop_pair = flux_prop_pair + 1e-9,
                             CN_dist = CN_dist + 1e-9), 
             x_var = "CN_dist", 
             y_var = "flux_prop_pair") +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5,
                                0.1),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1))
  

# Fig 3b
fig3b <- 
  basic_line(model = fluxprop_model_cnpdist, 
             dats = fluxes_long_perpair, 
             x_var = "CNP_dist", 
             y_var = "flux_prop_pair") +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01))

# Fig s1c
fig3c <- 
  basic_line(model = fluxraw_model_cndist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 0.0000000001), 
             x_var = "CN_dist", 
             y_var= "flux_raw_pair") +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5,
                                0.1),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_y_continuous(breaks = c(-22,
                                -20,
                                -10,
                                -5,
                                0,
                                5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1,
                                150)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))

# Fig s1d
fig3d <- 
  basic_line(model = fluxraw_model_cnpdist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 0.0000000001), 
             x_var = "CNP_dist", 
             y_var= "flux_raw_pair",
             scale = "log") +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01))+
  scale_y_continuous(breaks = c(-22,
                                -20,
                                -10,
                                -5,
                                0,
                                5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1,
                                150)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))

# Compile
fig3 <- 
  cowplot::plot_grid(plotlist = list(fig3a + 
                                       ggtitle("(a)"), 
                                     fig3b +
                                       ggtitle("(b)"),
                                     fig3c + 
                                       ggtitle("(c)"), 
                                     fig3d +
                                       ggtitle("(d)")),
                     nrow = 2)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig3.pdf"),
                fig3,
                height = 10,
                width = 10)




# Figure 4 (run script 03) ------------------------------------------------
# Fig 4a
fig4a <- 
  basic_line(model = totalflow_model_cndist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CN_dist", 
             y_var= "totalflow_Jy") +
  scale_x_continuous(breaks = c(-34,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1),
                     limits = c(-34, 0)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000))

# Fig 4b
fig4b <- 
  basic_line(model = predflow_model_cndist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CN_dist", 
             y_var= "predflow_Jy") +
  scale_x_continuous(breaks = c(-34,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1),
                     limits = c(-34, 0)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))
  

# Fig 4c
fig4c <- 
  basic_line(model = totalflow_model_cnpdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CNP_dist", 
             y_var = "totalflow_Jy") +
  scale_x_continuous(breaks = c(-34,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1),
                     limits = c(-34, 0)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000))

# Fig 4d
fig4d <- 
  basic_line(model = predflow_model_cnpdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CNP_dist", 
             y_var = "predflow_Jy") +
  scale_x_continuous(breaks = c(-34,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1),
                     limits = c(-34, 0)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))

# Compile
fig4 <- 
  cowplot::plot_grid(plotlist = list(fig4a +
                                       ggtitle("(a)"), 
                                     fig4b +
                                       ggtitle("(b)"), 
                                     fig4c +
                                       ggtitle("(c)"), 
                                     fig4d +
                                       ggtitle("(d)")),
                     nrow = 2)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig4.pdf"),
                fig4,
                height = 10,
                width = 10)



# Figure S1 (run script 03) --------
# Pair-scale but plotting with elemental ratios and functional groups 

# Prepare data
perpair_ratios <- 
  fluxes_long_perpair %>% 
  ## Get ratios
  dplyr::mutate(cons_CN = cons_C/cons_N,
                cons_CP = cons_C/cons_P,
                cons_NP = cons_N/cons_P,
                res_CN = res_C/res_N,
                res_CP = res_C/res_P,
                res_NP = res_N/res_P) %>% 
  ## Do comparisons
  dplyr::mutate(lower_cons_CN = ifelse(cons_CN < res_CN, "yes", "no"),
                lower_cons_CP = ifelse(cons_CP < res_CP, "yes", "no"),
                lower_cons_NP = ifelse(cons_NP < res_NP, "yes", "no")) 

# Flux proportion in CN plane
## Model effect
model_effect <- 
  brms::conditional_effects(fluxprop_model_cndist,
                            method = "fitted")[[1]] 
## Compare CN ratio
figs1a <- 
  ggplot(data = model_effect,
         aes(x = log(CN_dist),
             y = estimate__)) +
  ### Add raw data points
  geom_jitter(data = perpair_ratios %>% 
                dplyr::filter(!is.na(lower_cons_CN)),
              aes(x = log(CN_dist + 1e-9),
                  y = log(flux_prop_pair + 1e-9),
                  colour = lower_cons_CN),
              alpha = 0.4) +
  ### Colour of points
  scale_color_manual(values = c( "thistle", "thistle4"),
                     labels = c("No", "Yes"),
                     name = "Consumer CN \nlower than\nresource") +
  ### Line
  geom_line(linewidth = 1.2,
            colour = "steelblue") +
  ### Confidence interval
  geom_ribbon(aes(ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ### Styling
  xlab(full_label("CN_dist")) +
  ylab(full_label("flux_prop_pair")) + 
  ### Remove alpha from legend
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01))


## Compare CP ratios
figs1b <- 
  ggplot(data = model_effect,
         aes(x = log(CN_dist),
             y = estimate__)) +
  ### Add raw data points
  geom_jitter(data = perpair_ratios %>% 
                dplyr::filter(!is.na(lower_cons_CP)),
              aes(x = log(CN_dist + 1e-9),
                  y = log(flux_prop_pair + 1e-9),
                  colour = lower_cons_CP),
              alpha = 0.4) +
  ### Colour of points
  scale_color_manual(values = c( "thistle", "thistle4"),
                     labels = c("No", "Yes"),
                     name = "Consumer CP \nlower than\nresource") +
  ### Line
  geom_line(linewidth = 1.2,
            colour = "steelblue") +
  ### Confidence interval
  geom_ribbon(aes(ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ### Styling
  xlab(full_label("CN_dist")) +
  ylab(full_label("flux_prop_pair")) + 
  ### Remove alpha from legend
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01))

## Functional groups
figs1c <- 
  ggplot(data = model_effect,
         aes(x = log(CN_dist),
             y = estimate__)) +
  ### Add raw data points
  geom_jitter(data = perpair_ratios,
              aes(x = log(CN_dist + 1e-9),
                  y = log(flux_prop_pair + 1e-9),
                  colour = cons_functional_group),
              alpha = 0.4) +
  ### Colour of points
  scale_color_manual(values = c( "thistle4", "darkslategrey", "dodgerblue4", "lemonchiffon3","grey50"),
                     labels = c("Engulfer", "Gatherer", "Piercer", "Shredder", "Unknown"),
                     name = "Consumer \nfunctional group") +
  ### Line
  geom_line(linewidth = 1.2,
            colour = "steelblue") +
  ### Confidence interval
  geom_ribbon(aes(ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ### Styling
  xlab(full_label("CN_dist")) +
  ylab(full_label("flux_prop_pair")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ### Remove alpha from legend
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01))

# Raw fluxesin CN plane
## Model effect
model_effect <- 
  brms::conditional_effects(fluxraw_model_cndist,
                            method = "fitted")[[1]] 
## Compare CN ratio
figs1d <- 
  ggplot(data = model_effect,
         aes(x = log(CN_dist),
             y = estimate__)) +
  ### Add raw data points
  geom_jitter(data = perpair_ratios %>% 
                dplyr::filter(!is.na(lower_cons_CN)),
              aes(x = log(CN_dist + 1e-9),
                  y = log(flux_raw_pair + 1e-9),
                  colour = lower_cons_CN),
              alpha = 0.4) +
  ### Colour of points
  scale_color_manual(values = c( "thistle", "thistle4"),
                     labels = c("No", "Yes"),
                     name = "Consumer CN \nlower than\nresource") +
  ### Line
  geom_line(linewidth = 1.2,
            colour = "steelblue") +
  ### Confidence interval
  geom_ribbon(aes(ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ### Styling
  xlab(full_label("CN_dist")) +
  ylab(full_label("flux_raw_pair")) + 
  ### Remove alpha from legend
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))


## Compare CP ratios
figs1e <- 
  ggplot(data = model_effect,
         aes(x = log(CN_dist),
             y = estimate__)) +
  ### Add raw data points
  geom_jitter(data = perpair_ratios %>% 
                dplyr::filter(!is.na(lower_cons_CP)),
              aes(x = log(CN_dist + 1e-9),
                  y = log(flux_raw_pair + 1e-9),
                  colour = lower_cons_CP),
              alpha = 0.4) +
  ### Colour of points
  scale_color_manual(values = c( "thistle", "thistle4"),
                     labels = c("No", "Yes"),
                     name = "Consumer CP \nlower than\nresource") +
  ### Line
  geom_line(linewidth = 1.2,
            colour = "steelblue") +
  ### Confidence interval
  geom_ribbon(aes(ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ### Styling
  xlab(full_label("CN_dist")) +
  ylab(full_label("flux_raw_pair")) + 
  ### Remove alpha from legend
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01))

## Functional groups
figs1f <- 
  ggplot(data = model_effect,
         aes(x = log(CN_dist),
             y = estimate__)) +
  ### Add raw data points
  geom_jitter(data = perpair_ratios,
              aes(x = log(CN_dist + 1e-9),
                  y = log(flux_raw_pair + 1e-9),
                  colour = cons_functional_group),
              alpha = 0.4) +
  ### Colour of points
  scale_color_manual(values = c( "thistle4", "darkslategrey", "dodgerblue4", "lemonchiffon3","grey50"),
                     labels = c("Engulfer", "Gatherer", "Piercer", "Shredder", "Unknown"),
                     name = "Consumer \nfunctional group") +
  ### Line
  geom_line(linewidth = 1.2,
            colour = "steelblue") +
  ### Confidence interval
  geom_ribbon(aes(ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ### Styling
  xlab(full_label("CN_dist")) +
  ylab(full_label("flux_raw_pair")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ### Remove alpha from legend
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01,
                                1)) +
  scale_x_continuous(breaks = c(-32,
                                -20,
                                -10,
                                -5),
                     labels = c(0,
                                expression(paste("10"^{-10})),
                                expression(paste("10"^{-5})),
                                0.01))

# Combine plots
figs1 <- 
  cowplot::plot_grid(plotlist = list(figs1a +
                                       ggtitle("(a)"),
                                     figs1b +
                                       ggtitle("(b)"),
                                     figs1c +
                                       ggtitle("(c)"),
                                     figs1d +
                                       ggtitle("(d)"),
                                     figs1e +
                                       ggtitle("(e)"),
                                     figs1f +
                                       ggtitle("(f)")),
                     nrow = 2)
# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "figs1.pdf"),
                plot = figs1,
                width = 19,
                height = 10)



# Figures S2- S13 (run script 04)----------------------------------------------
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
  c("cons_CPI_crisp",
    "cons_P",
    "cons_CPI_crisp",
    "cons_P",
    "CN_dist",
    "CNP_dist",
    "CN_dist",
    "CNP_dist",
    "avg_CN_dist",
    "avg_CNP_dist",
    "avg_CN_dist",
    "avg_CNP_dist")
## Ylabs
ylabs <- 
  c("flux_prop",
    "flux_prop",
    "flux_raw",
    "flux_raw",
    "flux_prop_pair",
    "flux_prop_pair",
    "flux_raw_pair",
    "flux_raw_pair",
    "totalflow_Jy",
    "totalflow_Jy",
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
  if(i %in% c(1:4)){
    dats <- 
      fluxes_long_percons %>% 
      dplyr::filter(site == names(site_matrices)[i]) %>% 
      dplyr::mutate(flux_prop = flux_prop + 1e-3,
                    flux_raw = flux_raw + 1e-3)
  } else if(i %in% c(5:8)){
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
  
  
  ## Loop running through list of models
  for(j in 1:length(site_matrices)){
    ### Sanity check
    print(names(site_matrices)[j])
    
    ### Extract model
    model <- 
      sublist[[j]]
    
    ### Extract pd and percent in rope
    pd <- 
      bayestestR::pd(model)$pd[2]
    rope <- 
      bayestestR::rope(model)$ROPE_Percentage[2]
    
    ## Make check
    if(pd > 0.95 & rope < 0.05) 
      title <- paste(names(site_matrices)[j],
                     "***") else title <- names(site_matrices)[j]
    
    ## Make plot
    tempplot <- 
      basic_line(model = model, 
                 dats = model$data, 
                 x_var = xlabs[i], 
                 y_var = ylabs[i]) +
      ggtitle(title)
    
    # Add little space next to axis label where it gets cut off
    if(i %in% c(3, 4, 7, 8, 11, 12)){
      tempplot <- 
        tempplot +
        theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))
      
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
## fluxprop vs cpi
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[1]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs2.pdf"))
## fluxprop vs p
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[2]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs3.pdf")) 
## fluxraw vs cpi
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[3]] ,
                     nrow = 6,
                     ncol = 5) 
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs4.pdf")) 
## fluxraw vs p
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[4]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs5.pdf"))
## fluxprop vs cndist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[5]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs6.pdf"))
## fluxprop vs cnpdist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[6]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs7.pdf"))
## fluxraw vs cndist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[7]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs8.pdf"))
## fluxraw vs cnpdist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[8]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs9.pdf"))
## totalflow vs cn dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[9]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs10.pdf"))
## total flow vs cnpdist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[10]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs11.pdf"))
## pred flow vs cn dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[11]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs12.pdf"))
## pred flow vs cnp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[12]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 20,
                width = 20,
                filename = here::here("plots",
                                      "figs13.pdf"))





# Figure S14 (run script 04)--------------------------------------------------------------
# Fig 14a
figs14a <- 
  basic_line(model = totalflow_model_cndist_nozero, 
             dats = fluxes_long_perbrom %>% 
               dplyr::filter(avg_CN_dist > 0.001), 
             x_var = "avg_CN_dist", 
             y_var= "totalflow_Jy") +
  scale_x_continuous(breaks = c(-6.2,
                                -5,
                                -4,
                                -3,
                                0),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                1)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000))

# Fig 4b
figs14b <- 
  basic_line(model = predflow_model_cndist_nozero, 
             dats = fluxes_long_perbrom %>% 
               dplyr::filter(avg_CN_dist > 0.001), 
             x_var = "avg_CN_dist", 
             y_var= "predflow_Jy") +
  scale_x_continuous(breaks = c(-6.2,
                                -5,
                                -4,
                                -3,
                                0),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                1)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))


# Fig 4c
figs14c <- 
  basic_line(model = totalflow_model_cnpdist_nozero, 
             dats = fluxes_long_perbrom %>% 
               dplyr::filter(avg_CNP_dist > 0.001), 
             x_var = "avg_CNP_dist", 
             y_var= "totalflow_Jy") +
  scale_x_continuous(breaks = c(-6.2,
                                -5,
                                -4,
                                -3,
                                0),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                1)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000))

# Fig 4d
figs14d <- 
  basic_line(model = predflow_model_cnpdist_nozero, 
             dats = fluxes_long_perbrom  %>% 
               dplyr::filter(avg_CNP_dist > 0.001), 
             x_var = "avg_CNP_dist", 
             y_var= "predflow_Jy") +
  scale_x_continuous(breaks = c(-6.2,
                                -5,
                                -4,
                                -3,
                                0),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                1)) +
  scale_y_continuous(breaks = c(-5, 
                                -2.5, 
                                0, 
                                2.5, 
                                5, 
                                7),
                     labels = c(0.001, 
                                0.01, 
                                1, 
                                10, 
                                150,
                                1000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))


# Compile
figs14 <- 
  cowplot::plot_grid(plotlist = list(figs14a +
                                       ggtitle("(a)"), 
                                     figs14b +
                                       ggtitle("(b)"), 
                                     figs14c +
                                       ggtitle("(c)"), 
                                     figs14d +
                                       ggtitle("(d)")),
                     nrow = 2)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "figs14.pdf"),
                figs14,
                height = 10,
                width = 10)

