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
                            method = "fitted")$cons_P %>% 
  dplyr::mutate(estimate__ = exp(estimate__),
                lower__ = exp(lower__),
                upper__ = exp(upper__))

# Fig2a
fig2a <- 
  ggplot(data = fluxes_long_percons,
         aes(x = jitter(cons_P),
             y = cons_CPI_crisp)) +
  ## Add raw data points
  geom_jitter(colour = "lightsteelblue4",
              alpha = 0.4) +
  ## Line
  geom_line(data = model_effect,
            aes(x = cons_P,
                y = estimate__),
            linewidth = 1.2,
            colour = "steelblue") +
  ## Confidence interval
  geom_ribbon(data = model_effect,
              aes(x = cons_P,
                  ymin = lower__,
                  ymax = upper__),
              alpha = 0.5,
              colour = NA,
              fill = "steelblue") +
  ## Styling
  xlab(full_label("cons_P")) +
  ylab(full_label("cons_CPI_crisp")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

# Fig 2b
fig2b <- 
  basic_line(model = fluxprop_model_cpi, 
             dats = fluxes_long_percons, 
             x_var = "cons_CPI_crisp", 
             y_var= "flux_prop",
             scale = "sqrt")

# Fig 2c
fig2c <- 
  basic_line(model = fluxprop_model_p, 
             dats = fluxes_long_percons, 
             x_var = "cons_P", 
             y_var = "flux_prop",
             scale = "sqrt")

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
fig3a <- 
  basic_line(model = fluxprop_model_cndist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_prop_pair = flux_prop_pair + 0.0000000001), 
             x_var = "CN_dist", 
             y_var = "flux_prop_pair",
             scale = "log") +
  scale_y_continuous(trans = "log10",
                     breaks = c(0.0000000001,
                                0.00000001,
                                0.00001,
                                0.001,
                                0.1,
                                1),
                     labels = c("1 e-10",
                                "1 e-7",
                                "1 e-5",
                                "0.001",
                                "0.1",
                                "1"))

# Fig 3b
fig3b <- 
  basic_line(model = fluxprop_model_cnpdist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_prop_pair = flux_prop_pair + 0.0000000001), 
             x_var = "CNP_dist", 
             y_var = "flux_prop_pair",
             scale = "log") +
  scale_y_continuous(trans = "log10",
                     breaks = c(0.0000000001,
                                0.00000001,
                                0.00001,
                                0.001,
                                0.1,
                                1),
                     labels = c("1 e-10",
                                "1 e-7",
                                "1 e-5",
                                "0.001",
                                "0.1",
                                "1"))


# Compile
fig3 <- 
  cowplot::plot_grid(plotlist = list(fig3a + 
                                       ggtitle("(a)"), 
                                     fig3b +
                                       ggtitle("(b)")),
                     nrow = 1)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "fig3.pdf"),
                fig3,
                height = 5,
                width = 9)




# Figure 4 (run script 03) ------------------------------------------------
# Fig 4a
fig4a <- 
  basic_line(model = totalflow_model_CNdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CN_dist", 
             y_var= "totalflow_Jy",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 500, 1500),
                     labels = c("0.1", "1", "20", "500", "1500"))

# Fig 4b
fig4b <- 
  basic_line(model = predflow_model_CNdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CN_dist", 
             y_var= "predflow_Jy",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 500, 1500),
                     labels = c("0.1", "1", "20", "500", "1500"))

# Fig 4c
fig4c <- 
  basic_line(model = totalflow_model_CNPdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CNP_dist", 
             y_var = "totalflow_Jy",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 500, 1500),
                     labels = c("0.1", "1", "20", "500", "1500"))

# Fig 4d
fig4d <- 
  basic_line(model = predflow_model_CNPdist, 
             dats = fluxes_long_perbrom, 
             x_var = "avg_CNP_dist", 
             y_var = "predflow_Jy",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 500, 1500),
                     labels = c("0.1", "1", "20", "500", "1500"))

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



# Figure S1 (run script 03) ------------------------------------------------
# Fig s1a
figs1a <- 
  basic_line(model = fluxraw_model_cpi, 
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 0.0000000001), 
             x_var = "cons_CPI_crisp", 
             y_var= "flux_raw",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 300, 900),
                     labels = c("0.1", "1", "20", "300", "900")) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))

# Fig s1b
figs1b <- 
  basic_line(model = fluxraw_model_p, 
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 0.0000000001), 
             x_var = "cons_P", 
             y_var = "flux_raw",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 300, 900),
                     labels = c("0.1", "1", "20", "300", "900")) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))

# Fig s1c
figs1c <- 
  basic_line(model = fluxraw_model_cndist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 0.0000000001), 
             x_var = "CN_dist", 
             y_var= "flux_raw_pair",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 300, 900),
                     labels = c("0.1", "1", "20", "300", "900")) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))

# Fig s1d
figs1d <- 
  basic_line(model = fluxraw_model_cnpdist, 
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 0.0000000001), 
             x_var = "CNP_dist", 
             y_var= "flux_raw_pair",
             scale = "log") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, 1, 20, 300, 900),
                     labels = c("0.1", "1", "20", "300", "900")) +
  theme(axis.title.y = element_text(margin = margin(r = 0,  l = 5)))


# Compile
figs1 <- 
  cowplot::plot_grid(plotlist = list(figs1a +
                                       ggtitle("(a)"), 
                                     figs1b +
                                       ggtitle("(b)"),
                                     figs1c +
                                       ggtitle("(c)"),
                                     figs1d +
                                       ggtitle("(d)")), 
                     nrow = 2)

# Save
ggplot2::ggsave(filename = here::here("plots",
                                      "figs1.pdf"),
                figs1,
                height = 10,
                width = 12)


