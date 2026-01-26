# Get CPI and biomass
library(tidyverse)
library(here)
source(here::here("R_code",
                  "functions.R"))
library(hellometry)


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
# Olivier's paper on population in French Guiana, the complete cohort lasted 4 months
## Dézerald, O., C. Leroy, B. Corbara, A. Dejean, S. Talaga, and R. Céréghino. 
## 2017. Environmental drivers of invertebrate population dynamics in Neotropical 
## tank bromeliads. Freshw Biol 62: 229–242.
120

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
# Olivier's paper on Cyphon in French Guiana, the complete cohort lasted 6 months
## Dézerald, O., C. Leroy, B. Corbara, A. Dejean, S. Talaga, and R. Céréghino. 
## 2017. Environmental drivers of invertebrate population dynamics in Neotropical 
## tank bromeliads. Freshw Biol 62: 229–242.
180

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
                 cp3_val = c(115.5, 120, 333, 60.5, 
                             60, 369.1, 730, 450,
                             270, 200, 144, 180,
                             150, 
                             (60 + 676.5 + 730 + 450 + 270 + 200 + 180) / 7,
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
  ## Make defuzzified CPI column
  dplyr::rowwise() %>% 
  dplyr::mutate(CPI_defuzzified = ifelse(is.na(CP1)|is.na(CP2)|is.na(CP3),
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
                                "Hirudidae", family)) %>% 
  ## Fix functional groups that are NAs
  dplyr::mutate(functional_group = ifelse(bwg_name %in% c("Diptera.422","Odonata.11"),
                                               "engulfer",
                                               ifelse(bwg_name %in% c("Branchiopoda.4", "Diptera.536", "Diptera.537",
                                                                      "Diptera.440", "Diptera.441", "diptera.443"),
                                                      "filter.feeder",
                                                      ifelse(bwg_name %in% c("Diptera.419", "Diptera.405"),
                                                             "gatherer", 
                                                             ifelse(bwg_name == "Diptera.675",
                                                                    "shredder", functional_group))))) %>% 
  ## Finally make all NA consumers unknown functional group
  dplyr::mutate(functional_group = ifelse(is.na(functional_group),
                                               "uk", functional_group))


# Save
readr::write_csv(traits_CPI,
                 here::here("data",
                            "traits_CPI.csv"))


# Add biomass data --------------------------------------------------------
# Get abundance data from original database
biomass <- 
  fw_database %>% 
  purrr::pluck("abundance") %>% 
  ## Add taxonomy
  hellometry::add_taxonomy() %>% 
  ## Remove columns that are not that important yet
  dplyr::select(-dataset_id, -species_id, -domain, -kingdom,
                -phylum, -subphylum) %>% 
  ## We want dry estimates
  dplyr::mutate(biomass_type = "dry") %>% 
  ## Make column with unknown size class, size and that to receive biomass estimates
  dplyr::mutate(biomass_col = NA,
                size_col = "unknown",
                stage = ifelse(class == "Insecta",
                               "larva", "adult")) %>% 
  ## Add data from Dézerald et al. (2017, full ref above) to have estimates from Cyphon
  ## I used the mean weight from 10-20 individuals with the mean instar size 
  dplyr::bind_rows(tibble::tibble(bwg_name = rep("cypholivier", 10),
                                  class = rep("Insecta", 10),
                                  subclass = rep("Neoptera", 10),
                                  ord = rep("Coleoptera", 10),
                                  subord = rep("Polyphaga", 10),
                                  family = rep("Scirtidae", 10),
                                  genus = rep("Cyphon", 10),
                                  stage = rep("larva", 10),
                                  biomass_type = rep("dry", 10),
                                  biomass_col = c(0.00195, 0.0047, 0.0105, 0.0259, 
                                                  0.0598, 0.0975, 0.1644, 0.3386,
                                                  0.4649, 0.7296), ### NB denominator number of larvae in batch mass measurement
                                  size_col = c("0.14", "0.25", "0.34", "0.44", "0.54",
                                               "0.64", "0.74", "0.84", "0.94", "1.08"),
                                  abundance = c(451, 458, 557, 549, 448, 
                                                430, 338, 263, 221, 58))) %>% 
  ## Do the same for ostracods
  dplyr::bind_rows(tibble::tibble(bwg_name = rep("ostrolivier", 10),
                                  class = rep("Ostracoda", 10),
                                  subclass = rep("Podocopa", 10),
                                  ord = rep("Podocopida", 10),
                                  subord = rep("Cytherocopina", 10),
                                  family = rep("Limnocytheridae", 10),
                                  subfamily = rep("Timiriaseviinae", 10),
                                  genus = rep("Elpidium", 10),
                                  species = rep("Bromeliarum", 10),
                                  stage = rep("adult", 10),
                                  biomass_type = rep("dry", 10),
                                  biomass_col = c(0.0000543, 0.000096914, 0.000180958, 0.000335231, 
                                                  0.000714286, 0.0012, 0.002, 0.003363636,
                                                  0.0055, 0.0172), ### NB denominator number of larvae in batch mass measurement
                                  size_col = c("0.15", "0.24", "0.34", "0.44", "0.54",
                                               "0.64", "0.74", "0.84", "0.94", "1.02"),
                                  abundance = c(18, 141, 447, 765, 1010, 
                                                1779, 3189, 2108, 210, 8)))

# Suite of levels to be used for estimations
level_vec  <-
  c("bwg_name", "genus", "tribe", "subfamily", "family", "subord","ord", 
    "subclass", "class")

# Get the estimations
biomass <-
  hellometry(
    dats = biomass , 
    level_vec  = level_vec , 
    biomass_type = "dry",  
    use_BWG_db = TRUE, 
    no_BWG_data = FALSE,
    r_square_cutoff_upper = 1) ## model for ostracods has a lot of clear data

# Save
readr::write_csv(biomass$data,
                 here::here("data",
                            "dats_biomass.csv"))
