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
# Fig 2a
## Get effect
model_effect <- 
  brms::conditional_effects(cpi_model_P,
                            method = "fitted")$cons_P 
## Plot
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
## Plot
fig2b <- 
  basic_line(model = fluxraw_model_cpip, 
             index = 1,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-7), 
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

# Fig 2c
fig2c <- 
  basic_line(model = fluxraw_model_cpip,
             index = 2,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-7), 
             x_var = "cons_P", 
             y_var = "flux_raw") +
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
                                       ggtitle("(c)")),
                     nrow = 2)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig2.pdf"),
                fig2,
                height = 10,
                width = 10)

  



# Figure 3 (run script 03) ------------------------------------------------
# Fig 3a
## Main
fig3a1 <- 
  subset_line(model = fluxraw_model_cndist, 
             index = 3,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "CN_dist",
             y_var= "flux_raw_pair",
             colour = "lower_cons_CN") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 5)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 7)),
        axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.4)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 
## Inset
fig3a2 <- 
  basic_line(model = fluxraw_model_cndist, 
             index = 3,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "CN_dist",
             y_var= "flux_raw_pair",
             colour = "lower_cons_CN",
             scale = "none") +
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
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = -1,  
                             l = -1)) +
  geom_rect(xmin = -7,
            xmax = -1.9,
            ymin = -6.5,
            ymax = 3,
            fill = NA,
            colour = "black",
            linewidth = 1.5)
## Combine
fig3a <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(fig3a1) +
  cowplot::draw_plot(fig3a2, 
                     x = 0.11, 
                     y = .7, 
                     width = .3, 
                     height = .3)


# Fig 3b
## Main
fig3b1 <- 
  subset_line(model = fluxraw_model_cpdist, 
              index = 3,
              dats = fluxes_long_perpair %>% 
                dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
              x_var = "CP_dist",
              y_var= "flux_raw_pair",
              colour = "lower_cons_CP") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 5)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 7)),
        axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.4)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))  
## Inset
fig3b2 <- 
  basic_line(model = fluxraw_model_cpdist, 
              index = 3,
              dats = fluxes_long_perpair %>% 
                dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
              x_var = "CP_dist",
              y_var= "flux_raw_pair",
              colour = "lower_cons_CP",
              scale = "none") +
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
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0,  
                             l = 0)) +
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  geom_rect(xmin = -6,
            xmax = -1.9,
            ymin = -6.5,
            ymax = 3,
            fill = NA,
            colour = "black",
            linewidth = 1.5)
## Combine
fig3b <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(fig3b1) +
  cowplot::draw_plot(fig3b2, 
                     x = 0.11, 
                     y = .7, 
                     width = .3, 
                     height = .3)


# Fig 3c
## Main
fig3c1 <- 
  subset_line(model = fluxraw_model_cnpdist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "CNP_dist", 
             y_var= "flux_raw_pair") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 5)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 7)),
        axis.text = element_text(size = rel(1.4)),
        axis.title = element_text(size = rel(1.4))) 
## Inset
fig3c2 <- 
  basic_line(model = fluxraw_model_cnpdist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "CNP_dist", 
             y_var= "flux_raw_pair",
             scale = "none") +
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
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0,  
                             l = 0)) +
  xlab("") +
  ylab("")  +
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  geom_rect(xmin = -6.5,
            xmax = -1.9,
            ymin = -6.5,
            ymax = 3,
            fill = NA,
            colour = "grey30",
            linewidth = 1.5)
## Combine
fig3c <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(fig3c1) +
  cowplot::draw_plot(fig3c2, 
                     x = 0.11, 
                     y = .7, 
                     width = .3, 
                     height = .3) 

# Compile
fig3 <- 
  cowplot::plot_grid(plotlist = list(fig3a + 
                                       cowplot::draw_plot_label("(a)"), 
                                     fig3b + 
                                       cowplot::draw_plot_label("(b)"),
                                     fig3c + 
                                       cowplot::draw_plot_label("(c)")), 
                     nrow = 2)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig3.pdf"),
                fig3,
                height = 13,
                width = 20)




# Figure S1 (run script 03) ------------------------------------------------
# Fig S1a
## Main
figs1a1 <- 
  subset_line(model = totalflow_model_cndist, 
              index = 1,
              dats = fluxes_long_perbrom, 
              x_var = "avg_CN_dist",
              y_var= "totalflow_Jy") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 7)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 
## Inset
figs1a2 <- 
  basic_line(model = totalflow_model_cndist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CN_dist", 
             y_var= "totalflow_Jy",
             scale = "none") +
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
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = -1,  
                             l = -1)) +
  geom_rect(xmin = -7,
            xmax = -1.9,
            ymin = -2,
            ymax = 5.5,
            fill = NA,
            colour = "grey30",
            linewidth = 1.5)
## Combine
figs1a <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(figs1a1) +
  cowplot::draw_plot(figs1a2, 
                     x = 0.1, 
                     y = .7, 
                     width = .4, 
                     height = .3)

# Fig S1b
## Main
figs1b1 <- 
  subset_line(model = predflow_model_cndist, 
              index = 1,
              dats = fluxes_long_perbrom, 
              x_var = "avg_CN_dist",
              y_var= "predflow_Jy") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 8)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2))) 
## Inset
figs1b2 <- 
  basic_line(model = predflow_model_cndist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CN_dist", 
             y_var = "predflow_Jy",
             scale = "none") +
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
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = -1,  
                             l = -1)) +
  geom_rect(xmin = -7,
            xmax = -1.9,
            ymin = -5,
            ymax = 5,
            fill = NA,
            colour = "grey30",
            linewidth = 1.5)
## Combine
figs1b <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(figs1b1) +
  cowplot::draw_plot(figs1b2, 
                     x = 0.11, 
                     y = .7, 
                     width = .4, 
                     height = .3)

# Fig S1c
## Main
figs1c1 <- 
  subset_line(model = totalflow_model_cpdist, 
              index = 1,
              dats = fluxes_long_perbrom, 
              x_var = "avg_CP_dist",
              y_var= "totalflow_Jy") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 7)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 7)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 
## Inset
figs1c2 <- 
  basic_line(model = totalflow_model_cpdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CP_dist", 
             y_var= "totalflow_Jy",
             scale = "none") +
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
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = -1,  
                             l = -1)) +
  geom_rect(xmin = -7,
            xmax = -1.9,
            ymin = -2,
            ymax = 6,
            fill = NA,
            colour = "grey30",
            linewidth = 1.5)
## Combine
figs1c <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(figs1c1) +
  cowplot::draw_plot(figs1c2, 
                     x = 0.1, 
                     y = .7, 
                     width = .4, 
                     height = .3)

# Fig S1d
## Main
figs1d1 <- 
  subset_line(model = predflow_model_cpdist, 
              index = 1,
              dats = fluxes_long_perbrom, 
              x_var = "avg_CP_dist",
              y_var= "predflow_Jy") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 8)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2))) 
## Inset
figs1d2 <- 
  basic_line(model = predflow_model_cpdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CP_dist", 
             y_var = "predflow_Jy",
             scale = "none") +
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
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = -1,  
                             l = -1)) +
  geom_rect(xmin = -7,
            xmax = -1.9,
            ymin = -5,
            ymax = 5,
            fill = NA,
            colour = "grey30",
            linewidth = 1.5)
## Combine
figs1d <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(figs1d1) +
  cowplot::draw_plot(figs1d2, 
                     x = 0.11, 
                     y = .7, 
                     width = .4, 
                     height = .3)

# Fig S1e
## Main
figs1e1 <- 
  subset_line(model = totalflow_model_cnpdist, 
              index = 1,
              dats = fluxes_long_perbrom, 
              x_var = "avg_CNP_dist",
              y_var= "totalflow_Jy") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 7)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2))) 
## Inset
figs1e2 <- 
  basic_line(model = totalflow_model_cnpdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CNP_dist", 
             y_var= "totalflow_Jy",
             scale = "none") +
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
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = -1,  
                             l = -1)) +
  geom_rect(xmin = -7,
            xmax = -1.9,
            ymin = -2,
            ymax = 5.5,
            fill = NA,
            colour = "grey30",
            linewidth = 1.5)
## Combine
figs1e <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(figs1e1) +
  cowplot::draw_plot(figs1e2, 
                     x = 0.1, 
                     y = .7, 
                     width = .4, 
                     height = .3)

# Fig S1f
## Main
figs1f1 <- 
  subset_line(model = predflow_model_cnpdist, 
              index = 1,
              dats = fluxes_long_perbrom, 
              x_var = "avg_CNP_dist",
              y_var= "predflow_Jy") +
  scale_y_continuous(breaks = c(-5,
                                -2.5,
                                0,
                                2.5),
                     labels = c(0.005,
                                0.1,
                                1,
                                10),
                     limits = c(NA, 8)) +
  scale_x_continuous(breaks = c(seq(from = -6,
                                    to = -2,
                                    by = 1)),
                     labels = c(0.001,
                                0.005,
                                0.02,
                                0.05,
                                0.15)) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)),
        axis.title = element_text(size = rel(1.2))) 
## Inset
figs1f2 <- 
  basic_line(model = predflow_model_cnpdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CNP_dist", 
             y_var = "predflow_Jy",
             scale = "none") +
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
  scale_x_continuous(breaks = c(-30,
                                -20,
                                -10,
                                -5,
                                0),
                     labels = c(expression(paste("10"^{-15})),
                                expression(paste("10"^{-9})),
                                expression(paste("10"^{-5})),
                                expression(paste("10"^{-3})),
                                1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        plot.margin = margin(t = 0, 
                             r = 0, 
                             b = -1,  
                             l = -1)) +
  geom_rect(xmin = -7,
            xmax = -1.9,
            ymin = -5,
            ymax = 5,
            fill = NA,
            colour = "grey30",
            linewidth = 1.5)
## Combine
figs1f <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(figs1f1) +
  cowplot::draw_plot(figs1f2, 
                     x = 0.11, 
                     y = .7, 
                     width = .4, 
                     height = .3)
# Compile
figs1 <- 
  cowplot::plot_grid(plotlist = list(figs1a + 
                                       cowplot::draw_plot_label("(a)"), 
                                     figs1b +
                                       cowplot::draw_plot_label("(b)"), 
                                     figs1c +
                                       cowplot::draw_plot_label("(c)"), 
                                     figs1d +
                                       cowplot::draw_plot_label("(d)"), 
                                     figs1e +
                                       cowplot::draw_plot_label("(e)"), 
                                     figs1f +
                                       cowplot::draw_plot_label("(f)")),
                     nrow = 3)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "figs1.pdf"),
                figs1,
                height = 18,
                width = 15)



# Figures S2 - S12 (run script 04)----------------------------------------------
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

for(i in 4){

#for(i in 1:length(masterlist_sites)){
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
                                      "figs2.pdf")) 
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
                                      "figs3.pdf"))
## fluxraw vs cndist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[2]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs4.pdf"))
## fluxraw vs cpdist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[3]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs5.pdf"))
## fluxraw vs cnpdist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[1]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs6.pdf"))
## totalflow vs cn dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[5]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs7.pdf"))
## totalflow vs cp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[6]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs8.pdf"))
## totalflow vs cnp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[7]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs9.pdf"))
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
                                      "figs10.pdf"))
## pred flow vs cp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[9]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs11.pdf"))
## pred flow vs cnp dist
plot_grid <- 
  cowplot::plot_grid(plotlist = masterlist_siteplots[[10]],
                     nrow = 6,
                     ncol = 5)
ggplot2::ggsave(plot_grid,
                height = 15,
                width = 15,
                filename = here::here("plots",
                                      "figs12.pdf"))





