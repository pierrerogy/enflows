# Get CPI and stoichiometry
library(tidyverse)
library(bwgdata)
library(here)
library(fwdata)
source(here::here("R_code",
                  "functions.R"))
library(hellometry)
library(cowplot)
library(plotly)
library(brms)
library(bayestestR)

# Make rules for CPI ------------------------------------------------------
# Load fwdata
fw_database <- 
  readr::read_rds(here::here("data",
                             "scalewebs",
                             "fw_database.rds"))

# Extract organisms with non-0 values in CP3
## Problem here is that CP interval is not upper-bounded (>60 days)
## so we need to find estimate values for it
cp3 <- 
  fw_database$traits %>% 
  dplyr::select(bwg_name, dplyr::contains(c("ord",
                                            "subord",
                                            "family",
                                            "genus",
                                            "CP", ## CPI
                                            "FD", ## food
                                            "FG", ## feeding mode
                                            "functional_group"
  ))) %>% 
  ## Make a mock CPI column
  dplyr::filter(CP3  > 0)

# Literature search
# Trichoptera
# same genus but in Texas, median of collected larvae
## Norwood, J. C., and K. W. Stewart. 2002. Life history and 
## case-building behavior of Phylloicus ornatus (Trichoptera: Calamoceratidae) 
## in two spring-fed streams in Texas. Annals of the Entomological Society of 
## America 95: 44–56.
115.5

# Ostracoda
# order level, from three references on northern species, 6-24
# months, picked mean (15 months) and converted to days (* 30),
# value in upper range of life cycle of those in vernal ponds
## Smith, A. J., and L. D. Delorme. 2010. Ostracoda. In Ecology and 
## Classification of North American Freshwater Invertebrates. 
## pp. 725–771, Elsevier, Amsterdam, Netherlands.
450

# Odonata
# Leptagrion, "near one year", 333 days
## De Marco, P. Jr., and K. S. Furieri. 2000. Ecology of Leptagrion 
## perlongum Calvert, 1909: a bromeliad-dweller odonate species. 
## International Journal of Odonatology 11: 135–148.
# other references in Brazilian journals but no access to them
333

# Toxorhynchites
# many species in the genus, average of averages
## Steffan, W. A., and N. L. Evenhuis. 1981. Biology of Toxorhynchites. 
## Annu. Rev. Entomol. 26: 159–181.
# (41, 60.5, 21.4, 15.83, 36, 23.2, 21.5, 16.75) --> only value for CP3 is 60.5
60.5

# Curculionidae
# bromeliad-feeding weevil in Florida, 
# about two months to develop to adulthood, no precise measurement
## Cooper, T. M., and R. D. Cave. 2019. Metamasius callizona 
## (Coleoptera: Curculionidae): Fertility and larva survival to the 
## third instar in the laboratory. Florida Entomologist 102: 413-417.
60

# Dytiscidae
# on Copelatus, same genus than some bromeliad ones
# (369.1 ± 257.5 days, n = 32) for larvae to adult (307.4 ± 256.1 days, n =
# 32) for adult to death
## Watanabe, K., and S.-Y. Ohba. 2021. Physiological lifespan of 
## Copelatus parallelus Zimmermann, 1920 (Coleoptera: Dytiscidae)
## under laboratory conditions. The Coleopterists Bulletin 75: 516–518.
369.1

# Elateridae
# subfamilies in Mexico, including Agrypninae  of which some 
# species use bromeliads, then one brasilian species
## Zurita-García, M. L., P. Johnson, and S. Zaragoza-Caballero. 2014. 
## Biodiversidad de Elateridae (Coleoptera) en México. Revista Mexicana 
## de Biodiversidad 85: 303–311. --> 1-3 years as larvae, put 2
## Costa, C. 1970. Genus Pyrophorus 3. Life-history, larva and pupa 
## of Pyrophorus punctatissimus Blanchard (Col., Elateridae). Papéis 
## Avulsos de Zoologia 23: 69–76. --> 2 years as larva
730

# Elmidae
# all members of family, 6-24 months depending on species, take mean of 15
# mostly North American species
## Brown, H. P. 1987. Biology of riffle beetles. Annual 
## Review of Entomology 32: 253–273.
450

# Hydrophilidae
# southern brazil,same family than in database
# larvae predaceous (7-11 months to emergence, pick 9), adult detritivores (14 months)
## Albertoni, F. F., J. Steiner, and A. Zillikens. 2016. 
##The associated beetle fauna of Hohenbergia augusta and 
## Vriesea friburgensis (Bromeliaceae) in southern Brazil. 
## Journal of Natural History 50: 2917–2939.
270

# Lampyridae
# waiting reply from lead author
## Vaz, S., L. F. L. Da Silveira, and S. P. Rosa. 2020. Morphology and 
## life cycle of a new species of Psilocladus Blanchard, 1846 
## (Coleoptera, Lampyridae, Psilocladinae), the first known 
## bromeliad-inhabiting firefly. Pap. Avulsos Zool. 60: e202060(s.i.).24.
# so far last instar to emergence is 5 months, add 50 days after communication wiht author
200

# Planorbidae
# 9 growth curves in Guadeloupe, taking average of all 9, 20.5 weeks
## Leveque, C., and J. P. Pointier. 1976. Study of the growth of 
## Biomphalaria glabrata (Say) and other Planorbidae in Guadeloupe 
## (West Indies). Annals of Tropical Medicine & Parasitology 70: 199–204.
144

# Scirtidae
# some scirtes from Wisconsin raised in the laboratory, 11 months on average
## Kraatz, W. C. 1918. Scirtes tibialis, Guer., (Coleoptera-Dascyllidæ), 
## with observations on its life history. Annals of the Entomological 
## Society of America 11: 393–400.
# this one followed a population in the field, 11 months here too
## Lombardi, D. 1928. Contributo alla conoscenza dello Scirtes 
## hemisphaericus L. (Coleoptera-Helodidae). Bollettino del 
## Laboratorio di Entomologia del R. Istituto Superiore Agrario di 
## Bologna 1: 236–258.
333

# Tabanidae
# Several species from Minnesota, five months
## Philip, C. B. 1928. Methods of collecting and rearing the immature 
## stages of Tabanidae (Diptera). The Journal of Parasitology 14: 243–253.
150

# Now combine these values in CP3
## For order-level do averages of values in the group
cp3_search <- 
  tibble::tibble(taxon = c("Trichoptera", "Ostracoda", "Odonata", "Toxorhynchites", 
                           "Curculionidae", "Dytiscidae", "Elateridae", "Elmidae",
                           "Hydrophilidae", "Lampyridae", "Planorbidae", "Scirtidae",
                           "Tabanidae", "Coleoptera", "Diptera"),
                 cp3_val = c(115.5, 450, 333, 60.5, 
                             60, 369.1, 730, 450,
                             270, 200, 144, 333,
                             150, 
                             (60 + 676.5 + 730 + 450 + 270 + 200 + 333) / 7,
                             (60.5 + 150) / 2))


# Get CPI
traits_CPI <- 
  fw_database$traits %>% 
  dplyr::select(bwg_name, dplyr::contains(c("ord",
                                            "subord",
                                            "family",
                                            "genus",
                                            "CP", ## CPI
                                            "FD", ## food
                                            "FG", ## feeding mode
                                            "functional_group"))) %>% 
  ## Make crisp CPI column
  dplyr::rowwise() %>% 
  dplyr::mutate(CPI_crisp = ifelse(is.na(CP1)|is.na(CP2)|is.na(CP3),
                                   NA, mean(c(rep(10, CP1), rep(40, CP2), 
                                              rep(match_CP3(ord = ord, subord = subord, 
                                                            family = family, genus = genus,
                                                            cp3_search = cp3_search), 
                                                  CP3)))))  %>% 
  ## Ungroup
  dplyr::ungroup() %>% 
  ## Add some things and remove others
  dplyr::filter(!stringr::str_detect(string = bwg_name,
                                     pattern = "Acari|Nematoda")) %>% 
  dplyr::mutate(family = ifelse(stringr::str_detect(string = bwg_name,
                                                    pattern = "Hirud"),
                                "Hirudidae", family))


# Save
readr::write_csv(traits_CPI,
                 here::here("data",
                            "traits_CPI.csv"))




# Stoichiometry data -------------------------------------------------------------
# Data from 2017 paper
## https://www.frontiersin.org/articles/10.3389/fevo.2017.00110/full
## Raw data here obtained after communicating with lead author
## Data at the link above includes interpolated rows
stoich1 <- 
  readr::read_delim(here::here("data",
                               "Atacama.Cardoso.CR.data.csv")) %>% 
  ## Filter aquatic communities
  dplyr::filter(Habitat == "aquatic") %>% 
  ## Remove line column
  dplyr::select(-Lat, -Long, -Site, -Habitat, -CN, -CP, -NP, -TL, -id) %>% 
  ## Convert C and N percentage to proportion
  dplyr::mutate(C = C/100,
                N = N/100,
                P = P/100) %>% 
  ## Remove green species without match anywhere
  dplyr::filter(Genus.Species != "green") %>% 
  # Rename column to make it easier
  dplyr::rename(Sp = Genus.Species) %>% 
  ## Fix that typo in species name
  dplyr::mutate(Sp = ifelse(Sp == "Mecystogaster.modesta",
                                       "Mecistogaster.modesta", Sp)) %>% 
  ## Mutate Family column 
  dplyr::mutate(Family = ifelse(Sp == "Hirudinidae.sp",
                                "Hirudidae", Family)) %>% 
  ## Mutate Class column
  dplyr::mutate(Class = ifelse(Family == "Limnocytheridae",
                               "Ostracoda", Class)) %>%
  ## Mutate order column
  dplyr::mutate(Order = ifelse(Family == "Limnocytheridae",
                               "Podocopidae", Order)) %>% 
  ## Mutate Genus column 
  dplyr::mutate(Genus = ifelse(Sp == "Mecistogaster.modesta",
                               "Mecistogaster", Genus))

# Data from 2018 paper
## https://datadryad.org/stash/dataset/doi:10.5061/dryad.3cd2611
stoich2 <- 
  readr::read_csv(here::here("data",
                             "stoich2018.csv")) %>% 
  ## Remove extra columns
  dplyr::select(-TL, -CNmolar, -Ndet, -CNdet, -Tmin, - Tmax, -Tvar, -FG) %>% 
  ## Convert C and N percentage to proportion
  dplyr::mutate(C = C/100,
                N = N/100) %>% 
  ## Rename column
  dplyr::rename(BS.mg = Wmg) %>% 
  ## Fix typos
  dplyr::mutate(dplyr::across(dplyr::where(is.character),
                              ~stringr::str_replace(string = .,
                                                    pattern = "Ceratoponidae",
                                                    replacement = "Ceratopogonidae"))) %>% 
  dplyr::mutate(dplyr::across(dplyr::where(is.character),
                              ~stringr::str_replace(string = .,
                                                    pattern = "Dysticidae",
                                                    replacement = "Dytiscidae"))) %>% 
  dplyr::mutate(dplyr::across(dplyr::where(is.character),
                              ~stringr::str_replace(string = .,
                                                    pattern = "Forcipomya",
                                                    replacement = "Forcipomyia"))) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.character),
                              ~stringr::str_replace(string = .,
                                                    pattern = "Wyeomia",
                                                    replacement = "Wyeomyia"))) %>% 
  dplyr::mutate(dplyr::across(dplyr::where(is.character),
                              ~stringr::str_replace(string = .,
                                                    pattern = "Dolychopodidae",
                                                    replacement = "Dolichopodidae"))) %>% 
  dplyr::mutate(dplyr::across(dplyr::where(is.character),
                              ~stringr::str_replace(string = .,
                                                    pattern = "Brachicera",
                                                    replacement = "Brachycera")))


# Combine data
stoich_full <- 
  stoich1 %>% 
  ## Detect and remove duplicate rows in first tibble
  dplyr::anti_join(stoich1 %>% 
                     dplyr::select(Family, Sp, N, C) %>% 
                     dplyr::inner_join(stoich2 %>% 
                                         dplyr::select(Family, sp2 = Sp, N, C),
                                       by = c("Family", "N", "C")) %>% 
                     dplyr::select(-Family, -sp2),
                   by = c("Sp", "N", "C")) %>% 
  ## Adapt column name
  dplyr::rename(Site = Locality) %>% 
  ## Add second tibble
  dplyr::bind_rows(stoich2) %>% 
  ## Just do a quick reordering of columns
  dplyr::select(Site, Class, Order, Family, Genus, Sp, BS.mg, C, N, P) %>% 
  ## Give actual names needed by hellometry
  dplyr::rename(class = Class,
                ord = Order,
                family = Family,
                genus = Genus,
                species = Sp) %>% 
  ## Replace different character to later split whole column
  dplyr::mutate(species = stringr::str_replace(string = species,
                                              pattern = "\\.",
                                              replacement = "_")) %>% 
  ## Split species column
  tidyr::separate(col = species,
                  into = c("gen", "species"),
                  sep = "_") %>% 
  ## Remove extra splitted column
  dplyr::select(-gen) %>% 
  ## Fix typos
  dplyr::mutate(species = ifelse(species %in% c("bromeliarium", "bromelarium"),
                                 "bromeliarum",
                                 ifelse(species == "mel",
                                        "melanopus", species)),
                genus = ifelse(genus == "Elphidium",
                               "Elpidium", genus)) %>% 
  ## Add columns I need for hellometry
  dplyr::mutate(stage = ifelse(class == "Insecta",
                               "larva", "adult"),
                abundance = 1,
                bwg_name = ifelse(species == "davisi",
                                  "Diptera.113",
                                  ifelse(species == "andromache",
                                         "Odonata.5",
                                         ifelse(species == "elongatum",
                                                "Odonata.7",
                                                ifelse(genus == "Phylloicus",
                                                       "Trichoptera.1",
                                                       ifelse(genus == "Elpidium",
                                                              "Ostracoda.8",
                                                              ifelse(species == "kompii",
                                                                     "Diptera.11",
                                                                     ifelse(species == "neivai",
                                                                            "Diptera.12",
                                                                            ifelse(species == "jenningsi",
                                                                                   "Diptera.97",
                                                                                   ifelse(species == "rejector",
                                                                                          "Diptera.98",
                                                                                          ifelse(species == "modesta",
                                                                                                 "Odonata.4",
                                                                                                 ifelse(species == "melanopus",
                                                                                                        "Diptera.298",
                                                                                                        ifelse(species == "superterrenus",
                                                                                                               "Oligochaeta.13", NA))))))))))))) %>%
  ## Add the taxonomy columns
  hellometry::add_taxonomy() %>% 
  ## Tipulidae and Limoniidae are the same family now
  dplyr::mutate(family = ifelse(family %in% c("Tipulidae", "Limoniidae"),
                                "Tipulidae_Limoniidae", family))





# Explore a bit CN stoichiometric patterns with heatplots -----------------------------------
# Which taxa do we have
stoich_full %>% 
  dplyr::select(family) %>% 
  dplyr::group_by(family) %>% 
  dplyr::tally() %>% 
  print(n = 50)
stoich_full %>% 
  dplyr::select(genus) %>% 
  dplyr::group_by(genus) %>% 
  dplyr::tally() %>% 
  print(n = 50)
stoich_full %>% 
  dplyr::select(bwg_name) %>% 
  dplyr::group_by(bwg_name) %>% 
  dplyr::tally() %>% 
  print(n = 50)

# How many points per taxon
stoich_full %>% 
  dplyr::select(family) %>% 
  dplyr::group_by(family) %>% 
  dplyr::tally() %>% 
  print(n = 50)

# What is the spread of points in a taxa
## For each family
### Empty list
plotlist <- 
  list()
### Do plots for each family
for(i in unique(stoich_full$family)){
  plotlist[[length(plotlist) + 1]] <- 
    ggplot(stoich_full %>% 
             dplyr::filter(family == i), 
           aes(x = C, 
               y = N) ) +
    ylim(0.02, 0.16) +
    xlim(0.2, 0.8) +
    geom_hex(bins = 50) +
    ggtitle(i) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()
}
### Print plots together
familyplots <- 
  cowplot::plot_grid(plotlist = plotlist,
                     nrow = 4,
                     ncol = 6)
### Save plots
ggsave(here::here("plots",
                  "family_stoich.pdf"),
       familyplots,
       width = 25,
       height = 15)
## For each genus
### Empty list
plotlist <- 
  list()
### Do plots for each genus
for(i in unique(stoich_full$genus)){
  plotlist[[length(plotlist) + 1]] <- 
    ggplot(stoich_full %>% 
             dplyr::filter(genus == i), 
           aes(x = C, 
               y = N) ) +
    ylim(0.02, 0.16) +
    xlim(0.2, 0.8) +
    geom_hex(bins = 70) +
    ggtitle(i) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()
}
### Print plots together
genusplots <- 
  cowplot::plot_grid(plotlist = plotlist[2:length(plotlist)],
                     nrow = 4,
                     ncol = 7)
### Save plots
ggsave(here::here("plots",
                  "genus_stoich.pdf"),
       genusplots,
       width = 25,
       height = 15)
## For each bwg name
### Empty list
plotlist <- 
  list()
### Do plots for each genus
for(i in unique(stoich_full$bwg_name)){
  plotlist[[length(plotlist) + 1]] <- 
    ggplot(stoich_full %>% 
             dplyr::filter(bwg_name == i), 
           aes(x = C, 
               y = N) ) +
    ylim(0.02, 0.16) +
    xlim(0.2, 0.8) +
    geom_hex(bins = 70) +
    ggtitle(i) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()
}
### Print plots together
bwgnameplots <- 
  cowplot::plot_grid(plotlist = plotlist[2:length(plotlist)],
                     nrow = 4,
                     ncol = 3)
### Save plots
ggsave(here::here("plots",
                  "bwgname_stoich.pdf"),
       bwgnameplots,
       width = 10,
       height = 15)

# Explore a bit CN stoichiometric patterns with hulls -----------------------------------
# Empty plot
base <- 
  ggplot(stoich_full, 
         aes(x = C, 
             y = N))

# Make list of list of families to plot
taxon_list <- 
  list(c("Limnocytheridae", "Hirudidae",
         "Naididae", "Cytheridae"),
       c("Chironomidae", "Culicidae", "Tipulidae_Limoniidae",
         "Calamoceratidae", "Scirtidae", "Psychodidae",
         "Ceratopogonidae", "Syrphidae", "Sciaridae"),
       c("Corethrellidae", "Dytiscidae", "Coenagrionidae", 
         "Dolichopodidae", "Hydrophilidae", "Tabanidae",
         "Anisoptera", "Veliidae", "Pseudostigmatidae"))

# Loop to get family plots
for(i in 1:length(taxon_list)){
  ## Get hull for species in list
  hull <- 
    stoich_full %>%
    dplyr::filter(family %in% taxon_list[[i]]) %>% 
    dplyr::group_by(family) %>% 
    dplyr::select(family, C, N) %>% 
    dplyr::filter(!is.na(family) & !is.na(C) & !is.na(N)) %>% 
    dplyr::slice(grDevices::chull(C, N))
  
  ## Add hulls to plot
  assign(paste0("p", i),
         base + 
           aes(fill = factor(family)) +
           geom_polygon(data = hull, 
                     alpha = 0.5) +
           xlim(0.2, 0.75) +
           ylim(0.02, 0.16) +
           guides(fill = guide_legend(title = "Family")))
}

# Make list of list of genera to plot
taxon_list <- 
  list(c("Elpidium", "Aulophorus",
         "Dero"),
       c("Chironomus", "Culex", "Microculex",
         "Phylloicus", "Polypedilum", "Scirtes",
         "Anopheles", "Wyeomyia", "Copestilum",
         "Meromacrus", "Tentrepohlia", "Tanytarsus",
         "Telmatoscopus", "Cyphon", "Forcipomyia"),
       c("Corethrella", "Leptagrion", "Monopelopia", 
         "Mecistogaster", "Toxorhynchites", "Copelatus",
         "Fidena"))

# Loop to get genera plots
for(i in 1:length(taxon_list)){
  ## Get hull for species in list
  hull <- 
    stoich_full %>%
    dplyr::filter(genus %in% taxon_list[[i]]) %>% 
    dplyr::group_by(genus) %>% 
    dplyr::select(genus, C, N) %>% 
    dplyr::filter(!is.na(genus) & !is.na(C) & !is.na(N)) %>% 
    dplyr::slice(grDevices::chull(C, N))
  
  ## Add hulls to plot
  assign(paste0("p", i + 3),
         base + 
           aes(fill = factor(genus)) +
           geom_polygon(data = hull, 
                        alpha = 0.5) +
           xlim(0.2, 0.75) +
           ylim(0.02, 0.16) +
           guides(fill = guide_legend(title = "Genus")))
}

# Make list of list of species to plot
taxon_list <- 
  list(c("bromeliarum", "superterrenus"),
       c("davisi", "kompii", "neivai",
         "jenningsi", "rejector", "melanopus"),
       c("andromache", "elongatum", "modesta"))

# Loop to get species plots
for(i in 1:length(taxon_list)){
  ## Get hull for species in list
  hull <- 
    stoich_full %>%
    dplyr::filter(species %in% taxon_list[[i]]) %>% 
    dplyr::group_by(species) %>% 
    dplyr::select(species, C, N) %>% 
    dplyr::filter(!is.na(species) & !is.na(C) & !is.na(N)) %>% 
    dplyr::slice(grDevices::chull(C, N))
  
  ## Add hulls to plot
  assign(paste0("p", i + 6),
         base + 
           aes(fill = factor(species)) +
           geom_polygon(data = hull, 
                        alpha = 0.5) +
           xlim(0.2, 0.75) +
           ylim(0.02, 0.16) +
           guides(fill = guide_legend(title = "species")))
}



# Give title to plots in list
plotlist <- 
  list(p1 + 
         ggtitle("Non-insect families"),
       p2 + 
         ggtitle("Detritivore and mixed families"),
       p3 +
         ggtitle("Predatory families"),
       p4 + 
         ggtitle("Non-insects genera"),
       p5 + 
         ggtitle("Detritivore genera"),
       p6 +
         ggtitle("Predatory genera"),
       p7 + 
         ggtitle("Non-insect species"),
       p8 + 
         ggtitle("Detritivore species"),
       p9 +
         ggtitle("Predatory species"))

# Combine and save
## Combine
hullplots <- 
  cowplot::plot_grid(plotlist = plotlist,
                     nrow = 3,
                     ncol = 3)
## Save
ggsave(here::here("plots",
                  "hull_stoich.pdf"),
       hullplots,
       width = 15,
       height = 15)

# CNP scatterplot -----------------------------------
## Overall
plotly::save_image( 
  plot_ly(data = stoich_full,
        x = ~C, 
        y = ~N, 
        z = ~P, 
        type = "scatter3d", 
        mode = "markers",
        color = ~family) %>% 
  layout(legend = list(title = list(text = 'Family')),
         yaxis = list(range = list(0.02, 0.16)),
         xaxis = list(range = list(0.2, 0.8)),
         zaxis = list(range = list(0.001, 0.02))),
  file = here::here("plots", "plotly",
                    "all_families.png"))

## Do and save plots for each family
for(i in unique(stoich_full$family)){
  plotly::save_image( 
    plot_ly(data = stoich_full %>% 
              dplyr::filter(family == i),
            x = ~C, 
            y = ~N, 
            z = ~P, 
            type = "scatter3d", 
            mode = "markers") %>% 
      layout(yaxis = list(range = list(0.02, 0.16)),
             xaxis = list(range = list(0.2, 0.8)),
             zaxis = list(range = list(0.001, 0.02))),
    file = here::here("plots", "plotly",
                      paste0(i, ".png")))
    
}

## Do and save plots for each genus
for(i in unique(stoich_full$genus)){
  plotly::save_image( 
    plot_ly(data = stoich_full %>% 
              dplyr::filter(genus == i),
            x = ~C, 
            y = ~N, 
            z = ~P, 
            type = "scatter3d", 
            mode = "markers") %>% 
      layout(yaxis = list(range = list(0.02, 0.16)),
             xaxis = list(range = list(0.2, 0.8)),
             zaxis = list(range = list(0.001, 0.02))),
    file = here::here("plots", "plotly",
                      paste0(i, ".png")))
  
}

## Do and save plots for each species
for(i in unique(stoich_full$species)[!stringr::str_detect(string = unique(stoich_full$species),
                                                          pattern = "sp")]){
  plotly::save_image( 
    plot_ly(data = stoich_full %>% 
              dplyr::filter(species == i),
            x = ~C, 
            y = ~N, 
            z = ~P, 
            type = "scatter3d", 
            mode = "markers") %>% 
      layout(yaxis = list(range = list(0.02, 0.16)),
             xaxis = list(range = list(0.2, 0.8)),
             zaxis = list(range = list(0.001, 0.02))),
    file = here::here("plots", "plotly",
                      paste0(i, ".png")))
  
}


# Get models best explaining CNP for different taxa ---------------------
# Empty tibble to store models
model_pars <- 
  tibble::tibble()

# Run models
for(h in c("family", "genus", "bwg_name")){
  # Select list of taxa
  taxa <- 
    stoich_full %>%
    dplyr::select(taxa = dplyr::all_of(h)) %>%
    dplyr::filter(!is.na(taxa)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull()
  
  # Now for each taxa within this list
  for(i in taxa){
    ## Let's keep track of things
    print(h)
    print(i)
    ## For each element
    for(j in c("C", "N", "P")){
      ### Let's keep track of things
      print(j)
      ### Prepare data
      dats <- 
        stoich_full %>% 
        dplyr::filter(.data[[h]] == i) %>% 
        dplyr::select(Site, dplyr::all_of(h), element = dplyr::all_of(j), BS.mg) %>%
        dplyr::filter(!is.na(element) & !is.na(BS.mg))
      
      ### Count how many sites for random effect
      sites <- 
        length(unique(dats$Site))
      
      ### Count number of rows
      #### If less than 5 rows, add blank
      if(nrow(dats) < 3){
        model_pars <- 
          model_pars %>% 
          dplyr::bind_rows(tibble::tibble(taxon_level = h,
                                          taxon_name = i,
                                          element = j))
        
      } 
      #### If more than five rows, run model
      if(nrow(dats) > 3){
        ##### No random effects
        if(sites < 2){
          model <- 
            # brms::brm(element ~ log(BS.mg),
            #           iter = 2000,
            #           family = gaussian(),      
            #           control = list(adapt_delta = 0.99,
            #                          max_treedepth = 15),
            #           data = dats)
            lm(element ~ log(BS.mg),
               data = dats)
          
          }
        
        ##### Yes random effects
        if(sites >= 2){
          model <- 
            # brms::brm(element ~ (BS.mg) + (1|Site),
            #           iter = 2000,
            #           family = gaussian(),      
            #           control = list(adapt_delta = 0.99,
            #                          max_treedepth = 15),
            #           data = dats)
            lme4::lmer(element ~ log(BS.mg) + (1|Site),
                       data = dats)
          
          
          }
        
        ##### Extract model outputs
        outputs <- 
          # bayestestR::describe_posterior(model)
          car::Anova(model)
        
        summary <- 
          summary(model)$coefficients %>% 
          tibble::as_tibble()
        
        ##### Store in tibble
        model_pars <- 
          model_pars %>% 
          # dplyr::bind_rows(tibble::tibble(taxon_level = h,
          #                                 taxon_name = i,
          #                                 element = j,
          #                                 slope = outputs$Median[2],
          #                                 intercept = outputs$Median[1],
          #                                 pd = outputs$pd[2],
          #                                 prcnt_rope =outputs$ROPE_Percentage[2] ,
          #                                 r_squared = brms::bayes_R2(model)[1],
          #                                 nsites = sites,
          #                                 model = list(model)))
        dplyr::bind_rows(tibble::tibble(taxon_level = h,
                                        taxon_name = i,
                                        element = j,
                                        slope = summary$Estimate[2],
                                        intercept = summary$Estimate[1],
                                        p_value = outputs$`Pr(>Chisq)`,
                                        nsites = sites,
                                        model = list(model)))
      }
      
    }
    
  }

}

# Remove those models that did not work out
model_pars_fil <- 
  # model_pars %>% 
  # dplyr::filter(pd > 0.95 & prcnt_rope < 0.05) %>% 
  # dplyr::mutate(taxon = taxon_name)
  model_pars %>% 
  dplyr::filter(p_value < 0.05) %>% 
  dplyr::mutate(taxon = taxon_name)

# Use biomass data to estimate stoichiometry -----------------------------------------------------
## Load equation file
equations <- 
  readr::read_csv(here::here("data",
                             "equations.csv")) %>% 
  dplyr::select(element:slope)


## Get lowest taxonomic level for each species in database
last_taxon <-
  traits_CPI %>% 
  ## Tipulidae and Limoniidae are the same family now
  dplyr::mutate(family = ifelse(family %in% c("Tipulidae", "Limoniidae"),
                                "Tipulidae_Limoniidae", family)) %>% 
  ## Remove columns I do not want,
  dplyr::select(-dplyr::contains("FD"), 
                -dplyr::contains("FG")) %>% 
  ## Very annoyed here
  dplyr::mutate(genus = ifelse(is.na(genus),
                               family, genus)) %>% 
  ## Make column with last taxon available
  dplyr::mutate(last_taxon = ifelse(!is.na(genus),
                                     genus, ifelse(!is.na(family),
                                                   family, NA))) %>% 
  ## Make some exceptions
  dplyr::mutate(last_taxon = ifelse(subfamily == "Orthocladiinae",
                                    "Orthocladiinae",
                                    ifelse(family == "Chironomidae" &
                                             genus %notin% c("Chironomus", "Monopelopia",
                                                             "Chironomidae"), 
                                           "Chironomidae",
                                           ifelse(family == "Culicidae" &
                                                    genus %notin% c("Anopheles", "Culex", "Microculex",
                                                                    "Toxorhynchites", "Wyeomyia", NA),
                                                  "Culicidae", last_taxon)))) %>%
  ## Last attempt
  dplyr::mutate(last_taxon = ifelse(family %in% c("Syrphidae", "Tabanidae",
                                                   "Scirtidae", "Tipulidae_Limoniidae",
                                                   "Limnocytheridae",
                                                   "Muscidae", "Corethrellidae",
                                                   "Coenagrionidae", "Dytiscidae",
                                                   "Pseudostigmatidae", "Psychodidae"),
                                     family, last_taxon)) %>%
  ## Get names and taxa
  dplyr::select(bwg_name, taxon = last_taxon)

## Read biomass data
biomass_stoich  <- 
  readr::read_csv(here::here("data",
                           "biomass_clean.csv")) %>% 
  ## Add last taxon
  dplyr::left_join(last_taxon,
                   by = "bwg_name") %>% 
  ## Add C
  dplyr::left_join(equations %>% 
                     dplyr::filter(element == "C"),
                   by = "taxon") %>% 
  dplyr::mutate(C = intercept + slope *(avg_bodymass_g/1000)) %>% 
  dplyr::select(-slope, -element, -intercept) %>% 
  ## Add N
  dplyr::left_join(equations %>% 
                     dplyr::filter(element == "N"),
                   by = "taxon") %>% 
  dplyr::mutate(N = intercept + slope * (avg_bodymass_g/1000)) %>% 
  dplyr::select(-slope, -element, -intercept) %>% 
  ## Add P
  dplyr::left_join(equations %>% 
                     dplyr::filter(element == "P"),
                   by = "taxon") %>% 
  dplyr::mutate(P = intercept + slope * (avg_bodymass_g/1000)) %>% 
  dplyr::select(-slope, -element, -intercept) %>% 
  ## Remove 0s put on proportion scale
  dplyr::mutate(dplyr::across(c(C, N, P),
                              ~ ifelse(. < 0, 
                                       NA, ./100)))



## Save data
readr::write_csv(biomass_stoich,
                 here::here("data",
                            "biomass_stoich.csv"))
















# Try to find the same equations than Ang ---------------------------------
# Get models best explaining CNP for different taxa ---------------------
# Empty tibble to store models
model_pars <- 
  tibble::tibble()


test <- 
  stoich1 %>% 
  ## Adapt column name
  dplyr::rename(Site = Locality) %>% 
  ## Just do a quick reordering of columns
  dplyr::select(Site, Class, Order, Family, Genus, Sp, BS.mg, C, N, P) %>% 
  ## Give actual names needed by hellometry
  dplyr::rename(class = Class,
                ord = Order,
                family = Family,
                genus = Genus,
                species = Sp) %>% 
  ## Replace different character to later split whole column
  dplyr::mutate(species = stringr::str_replace(string = species,
                                               pattern = "\\.",
                                               replacement = "_")) %>% 
  ## Split species column
  tidyr::separate(col = species,
                  into = c("gen", "species"),
                  sep = "_") %>% 
  ## Remove extra splitted column
  dplyr::select(-gen) %>% 
  ## Fix typos
  dplyr::mutate(species = ifelse(species %in% c("bromeliarium", "bromelarium"),
                                 "bromeliarum",
                                 ifelse(species == "mel",
                                        "melanopus", species)),
                genus = ifelse(genus == "Elphidium",
                               "Elpidium", genus)) %>% 
  ## Add columns I need for hellometry
  dplyr::mutate(stage = ifelse(class == "Insecta",
                               "larva", "adult"),
                abundance = 1,
                bwg_name = ifelse(species == "davisi",
                                  "Diptera.113",
                                  ifelse(species == "andromache",
                                         "Odonata.5",
                                         ifelse(species == "elongatum",
                                                "Odonata.7",
                                                ifelse(genus == "Phylloicus",
                                                       "Trichoptera.1",
                                                       ifelse(genus == "Elpidium",
                                                              "Ostracoda.8",
                                                              ifelse(species == "kompii",
                                                                     "Diptera.11",
                                                                     ifelse(species == "neivai",
                                                                            "Diptera.12",
                                                                            ifelse(species == "jenningsi",
                                                                                   "Diptera.97",
                                                                                   ifelse(species == "rejector",
                                                                                          "Diptera.98",
                                                                                          ifelse(species == "modesta",
                                                                                                 "Odonata.4",
                                                                                                 ifelse(species == "melanopus",
                                                                                                        "Diptera.298",
                                                                                                        ifelse(species == "superterrenus",
                                                                                                               "Oligochaeta.13", NA))))))))))))) %>%
  ## Add the taxonomy columns
  hellometry::add_taxonomy() %>% 
  ## Tipulidae and Limoniidae are the same family now
  dplyr::mutate(family = ifelse(family %in% c("Tipulidae", "Limoniidae"),
                                "Tipulidae_Limoniidae", family))

# Run models
for(h in c("family", "genus", "bwg_name")){
  # Select list of taxa
  taxa <- 
    test %>%
    dplyr::select(taxa = dplyr::all_of(h)) %>%
    dplyr::filter(!is.na(taxa)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull()
  
  # Now for each taxa within this list
  for(i in taxa){
    ## Let's keep track of things
    print(h)
    print(i)
    ## For each element
    for(j in c("C", "N", "P")){
      ### Let's keep track of things
      print(j)
      ### Prepare data
      dats <- 
        test %>% 
        dplyr::filter(.data[[h]] == i) %>% 
        dplyr::select(Site, dplyr::all_of(h), element = dplyr::all_of(j), BS.mg) %>%
        dplyr::filter(!is.na(element) & !is.na(BS.mg))
      
      ### Count how many sites for random effect
      sites <- 
        length(unique(dats$Site))
      
      ### Count number of rows
      #### If less than 5 rows, add blank
      if(nrow(dats) < 3){
        model_pars <- 
          model_pars %>% 
          dplyr::bind_rows(tibble::tibble(taxon_level = h,
                                          taxon_name = i,
                                          element = j))
        
      } 
      #### If more than five rows, run model
      if(nrow(dats) > 3){
        ##### No random effects
        if(sites < 2){
          model <- 
            # brms::brm(log10(element) ~ log10(BS.mg),
            #           iter = 2000,
            #           family = gaussian(),      
            #           control = list(adapt_delta = 0.99,
            #                          max_treedepth = 15),
            #           data = dats)
            lm(log(element) ~ log(BS.mg),
               data = dats)
          
        }
        
        ##### Yes random effects
        if(sites >= 2){
          model <- 
            # brms::brm(log10(element) ~ log10(BS.mg) + (1|Site),
            #           iter = 2000,
            #           family = gaussian(),      
            #           control = list(adapt_delta = 0.99,
            #                          max_treedepth = 15),
            #           data = dats)
            lme4::lmer(log(element) ~ log(BS.mg) + (1|Site),
                       data = dats)
          
          
        }
        
        ##### Extract model outputs
        outputs <- 
          # bayestestR::describe_posterior(model)
          car::Anova(model)
        
        summary <- 
          summary(model)$coefficients %>% 
          tibble::as_tibble()
        
        ##### Store in tibble
        model_pars <- 
          model_pars %>% 
          # dplyr::bind_rows(tibble::tibble(taxon_level = h,
          #                                 taxon_name = i,
          #                                 element = j,
          #                                 slope = outputs$Median[2],
          #                                 intercept = outputs$Median[1],
          #                                 pd = outputs$pd[2],
          #                                 prcnt_rope =outputs$ROPE_Percentage[2] ,
          #                                 r_squared = brms::bayes_R2(model)[1],
          #                                 nsites = sites,
          #                                 model = list(model)))
          dplyr::bind_rows(tibble::tibble(taxon_level = h,
                                          taxon_name = i,
                                          element = j,
                                          slope = summary$Estimate[2],
                                          intercept = summary$Estimate[1],
                                          p_value = outputs$`Pr(>Chisq)`,
                                          nsites = sites,
                                          model = list(model)))
      }
      
    }
    
  }
  
}

# Remove those models that did not work out
model_pars_fil <- 
  # model_pars %>% 
  # dplyr::filter(pd > 0.95 & prcnt_rope < 0.05) %>% 
  # dplyr::mutate(taxon = taxon_name)
  model_pars %>% 
  dplyr::filter(p_value < 0.01) %>% 
  dplyr::mutate(taxon = taxon_name)



# Sort later  SAME THING WITH WHICH EQUATIONS TO USE--------------------------------------------------------------


## Add C
biomass_stoich <- 
  biomass_data %>% 
  ## Add C
  dplyr::rename(BS.mg = avg_bodymass_g) %>%
  do_predictions(.,
                 model_list = model_pars %>% 
                   dplyr::filter(element == "C")) %>% 
  dplyr::rename(C = Estimate,
                C_err = Est.Error,
                C_low = Q2.5,
                C_high = Q97.5)

## Add N
biomass_stoich2 <- 
  biomass_stoich %>% 
  ## Add N
  do_predictions(.,
                 model_list = model_pars_fil %>% 
                   dplyr::filter(element == "N")) %>% 
  dplyr::rename(N = Estimate,
                N_err = Est.Error,
                N_low = Q2.5,
                N_high = Q97.5)

## Add P
biomass_stoich3 <- 
  biomass_stoich2 %>% 
  ## Add N
  do_predictions(do_predictions(.,
                                model_list = model_pars_fil %>% 
                                  dplyr::filter(element == "P")) %>% 
  dplyr::rename(P = Estimate,
                P_err = Est.Error,
                P_low = Q2.5,
                P_high = Q97.5)
  

## Save

tt=
  biomass_data %>% 
  ## Add C
  dplyr::left_join(stoich_dic %>% 
                     dplyr::filter(element == "C") %>% 
                     dplyr::select(-taxon, -element),
                   by = "bwg_name") %>% 
  dplyr::rename(BS.mg = avg_bodymass_g) %>%
  do_predictions() %>% 
  dplyr::rename(C = Estimate,
                C_err = Est.Error,
                C_low = Q2.5,
                C_high = Q97.5)
