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

# Load binary site matrices
site_matrices <- 
  load_matrices(folder = here::here("data",
                                    "matrices"),
                ending = "_PM.txt") %>% 
  make_binary_matrices(matrix_list = .,
                       cutoff = 0.02)

# Load biomass data
biomass_clean <- 
  readr::read_csv(here::here("data",
                             "biomass_clean.csv"))



# Load models
## List of models
model_list <- 
  c(
    "consflow_model_axis1",
    "consflow_model_axis2",
    "consflow_model_axis3",
    "consflow_model_axis4",
    "pairflow_model_overall",
    "pairflow_model_axis1",
    "pairflow_model_axis2",
    "pairflow_model_axis3",
    "pairflow_model_axis4",
    "totalflow_model_vol",
    "gini_model_vol",
    "totalflow_model_dist",
    "gini_model_dist"
    )

##  Check if model present in environment, and if not load it
for (model in model_list) {
  if (!exists(model)) {
    assign(model, 
           readr::read_rds(here::here("data",
                                      paste0(model, ".RDS"))))}}




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
  site_matrices[[which(names(site_matrices) == "PTL")]]

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
  dplyr::left_join(biomass_clean %>% 
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
                           "pitiweb.pdf"),
                pitiweb,
                height = 10,
                width = 10)
ggplot2::ggsave(here::here("plots",
                           "pitiweb.jpeg"),
                pitiweb,
                height = 10,
                width = 10)

# Figure 2 (run script 05) ------------------------------------------------
# Fig 2a
fig2a <- 
  basic_line(model = consflow_model_axis1, 
             index = 1,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-3), 
             x_var = "Axis.1", 
             y_var= "flux_raw",
             scale = "x_scale") +
  scale_y_continuous(limits = c(-7, 10),
                     breaks = c(-7,
                                seq(from = -5, 
                                    to = 10,
                                    by = 5)),
                     labels= c(0,
                               0.001,
                               1,
                               150,
                               2000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))

# Fig 2b
fig2b <- 
  basic_line(model = consflow_model_axis2, 
             index = 1,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-3), 
             x_var = "Axis.2", 
             y_var= "flux_raw",
             scale = "x_scale") +
  scale_y_continuous(limits = c(-7, 10),
                     breaks = c(-7,
                                seq(from = -5, 
                                    to = 10,
                                    by = 5)),
                     labels= c(0,
                               0.001,
                               1,
                               150,
                               2000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))

# Fig 2c
fig2c <- 
  basic_line(model = consflow_model_axis3, 
             index = 1,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-3), 
             x_var = "Axis.3", 
             y_var= "flux_raw",
             scale = "x_scale") +
  scale_y_continuous(limits = c(-7, 10),
                     breaks = c(-7,
                                seq(from = -5, 
                                    to = 10,
                                    by = 5)),
                     labels= c(0,
                               0.001,
                               1,
                               150,
                               2000)) +
  theme(axis.title.y = element_text(margin = margin(r = 0, l = 5)))

# Fig 2d
fig2d <- 
  basic_line(model = consflow_model_axis4, 
             index = 1,
             dats = fluxes_long_percons %>% 
               dplyr::mutate(flux_raw = flux_raw + 1e-3), 
             x_var = "Axis.4", 
             y_var= "flux_raw",
             scale = "x_scale") +
  scale_y_continuous(limits = c(-7, 10),
                     breaks = c(-7,
                                seq(from = -5, 
                                    to = 10,
                                    by = 5)),
                     labels= c(0,
                               0.001,
                               1,
                               150,
                               2000)) +
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
                                       ggtitle("(d)")),
                     nrow = 2,
                     ncol = 2)

# Save
## PDF
ggplot2::ggsave(filename = here::here("plots",
                                      "fig2.pdf"),
                fig2,
                height = 8,
                width = 9)
## PNG
ggplot2::ggsave(filename = here::here("plots",
                                      "fig2.png"),
                fig2,
                bg = "white",
                height = 8,
                width = 9)

  



# Figure 3 (run script 05) ------------------------------------------------
# Fig 3a
fig3a <- 
  basic_line(model = pairflow_model_axis1, 
             index = 1,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "dist_axis1",
             y_var = "flux_raw_pair",
             scale = "x_scale") +
  scale_y_continuous(breaks = c(-7,
                                -5,
                                -2.5,
                                0,
                                2.5,
                                5),
                     labels = c(0,
                                0.005,
                                0.1,
                                1,
                                10, 150),
                     limits = c(NA, 5)) +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 10)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 3b
fig3b <- 
  basic_line(model = pairflow_model_axis2, 
             index = 1,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "dist_axis2",
             y_var = "flux_raw_pair",
             scale = "x_scale") +
  scale_y_continuous(breaks = c(-7,
                                -5,
                                -2.5,
                                0,
                                2.5,
                                5),
                     labels = c(0,
                                0.005,
                                0.1,
                                1,
                                10, 150),
                     limits = c(NA, 5)) +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 10)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 3c
fig3c <- 
  basic_line(model = pairflow_model_axis3, 
             index = 1,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "dist_axis3",
             y_var = "flux_raw_pair",
             scale = "x_scale") +
  scale_y_continuous(breaks = c(-7,
                                -5,
                                -2.5,
                                0,
                                2.5,
                                5),
                     labels = c(0,
                                0.005,
                                0.1,
                                1,
                                10, 150),
                     limits = c(NA, 5)) +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 10)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 3d
fig3d <- 
  basic_line(model = pairflow_model_axis4, 
             index = 1,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "dist_axis4",
             y_var = "flux_raw_pair",
             scale = "x_scale") +
  scale_y_continuous(breaks = c(-7,
                                -5,
                                -2.5,
                                0,
                                2.5,
                                5),
                     labels = c(0,
                                0.005,
                                0.1,
                                1,
                                10, 150),
                     limits = c(NA, 5)) +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 10)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 3e
fig3e <- 
  basic_line(model = pairflow_model_overall, 
             index = 1,
             dats = fluxes_long_perpair %>% 
               dplyr::mutate(flux_raw_pair = flux_raw_pair + 1e-3), 
             x_var = "func_dist",
             y_var = "flux_raw_pair",
             scale = "x_scale") +
  scale_y_continuous(breaks = c(-7,
                                -5,
                                -2.5,
                                0,
                                2.5,
                                5),
                     labels = c(0,
                                0.005,
                                0.1,
                                1,
                                10, 150),
                     limits = c(NA, 5)) +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 10)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Compile
fig3 <- 
  cowplot::plot_grid(plotlist = list(fig3a +
                                       ggtitle("(a)"), 
                                     fig3b +
                                       ggtitle("(b)"),
                                     fig3c +
                                       ggtitle("(c)"), 
                                     fig3d +
                                       ggtitle("(d)"),
                                     fig3e +
                                       ggtitle("(e)")),
                     nrow = 3,
                     ncol = 2)

# Save
## PDF
ggplot2::ggsave(filename = here::here("plots",
                                      "fig3.pdf"),
                fig3,
                height = 12,
                width = 9)
## PNG
ggplot2::ggsave(filename = here::here("plots",
                                      "fig3.png"),
                fig3,
                bg = "white",
                height = 12,
                width = 9)

# Figure 4 (run script 05)  ------------------------------------------------
# Fig 4a
fig4a <- 
  basic_line(model = totalflow_model_dist, 
             index = 1,
             dats = fluxes_long_perbrom, 
             x_var = "avg_func_dist",
             y_var= "totalflow_Jy",
             scale = "x_scale") +
  scale_y_continuous(breaks = c(-3,
                                0,
                                2.5,
                                5,
                                7),
                     labels = c(0,
                                1,
                                10,
                                100,
                                1000),
                     limits = c(-3.2, 10)) +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 5)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4b
fig4c <- 
  basic_line(model = gini_model_dist, 
             index = 1,
             dats = fluxes_long_perbrom, 
             x_var = "avg_func_dist",
             y_var= "gini",
             scale = "both_original") +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 5)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4c
fig4d <- 
  basic_line(model = totalflow_model_vol, 
             index = 1,
             dats = fluxes_long_perbrom, 
             x_var = "prop_hull",
             y_var= "totalflow_Jy",
             scale = "none") +
  scale_y_continuous(breaks = c(-3,
                                0,
                                2.5,
                                5,
                                7),
                     labels = c(0,
                                1,
                                10,
                                100,
                                1000),
                     limits = c(-3.2, 10)) +
  scale_x_continuous(breaks = c(seq(from = -15,
                                    to = 0,
                                    by = 3)),
                     labels = c(expression(paste("10"^{-7})),
                                expression(paste("10"^{-6})),
                                expression(paste("10"^{-4})),
                                0.002,
                                0.05,
                                1)) +
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 5)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) 

# Fig 4d
fig4d <- 
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
  theme(axis.title.y = element_text(margin = margin(r = 2,  l = 5)),
        axis.title.x = element_text(margin = margin(t = 15)),
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
                                       ggtitle("(d)")),
                     nrow = 2)

# Save
## PDF
ggplot2::ggsave(filename = here::here("plots",
                                      "fig4.pdf"),
                fig4,
                height = 8,
                width = 9)
## PNG
ggplot2::ggsave(filename = here::here("plots",
                                      "fig4.png"),
                fig4,
                bg = "white",
                height = 8,
                width = 9)


