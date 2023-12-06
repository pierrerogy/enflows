# Get raw FPOM CPOM
library(tidyverse)
library(bwgdata)
library(here)
library(bwgdata)
library(fwdata)
source(here::here("R_code",
                  "functions.R"))

# Access BWG detritus data ----------
## Custom function to access the bwg dataset, because url changed
credentials <- new.env(parent = emptyenv())
bwg_auth <- function() {
  ## should check if it exists, and not rerun the function.
  
  username <- readline("USERNAME: ")
  password <- getPass::getPass("PASSWORD: ")
  
  
  ## create POST 
  url <- "https://www.zoology.ubc.ca/~srivast/bwgdb/v2/?route=users&action=login#/dashboard"
  
  body <- list(username = username,
               password = password)
  
  r <- httr::POST(url, body = body, encode = "json")
  
  ## check status code
  httr::stop_for_status(r)
  
  answer <- httr::content(r, as = "text")
  
  credentials$token  <- jsonlite::fromJSON(answer)[["results"]]
  
}

## Get data
bwgdet <- 
  ## Get data
  bwgdata::bwg_get("bromeliads") %>% 
  ## Select detritus columns
  dplyr::select(bromeliad_id, visit_id, detritus) %>%
  ## Unnest detritus columns
  ### Here we will have a row for each detritus size range (min, max and mass)
  tidyr::unnest(detritus) %>% 
  ## Add that mass is in g (based on the website)
  dplyr::rename(mass_g = mass) %>%  
  ## Make numeric
  dplyr::mutate(dplyr::across(c(min, max, mass_g),
                              ~as.numeric(.))) %>% 
  ## Assign to categories
  ### First replace NAs with arbitrary value absent from original ranges
  dplyr::mutate(max = ifelse(is.na(max),
                             4000, max)) %>% 
  ### Pick 1000 as cutoff between CPOM and FPOM, to maximise assignment to a category
  dplyr::mutate(cat = ifelse(min >= 1000,
                             "coarse", ifelse(max <= 1000,
                                              "fine", NA))) %>% 
  ### Some with group all others under all
  dplyr::mutate(cat = ifelse(is.na(cat),
                             "all", cat)) %>% 
  ## Remove max and moin
  ## Now sum up categories for each bromeliad
  ## Keeping the lowest minimum and the highest maximum
  dplyr::group_by(bromeliad_id, visit_id, cat) %>% 
  dplyr::summarise(dplyr::across(mass_g, sum),
                   dplyr::across(min, min),
                   dplyr::across(max, max)) %>% 
  dplyr::ungroup()

# Look how detritus is split across all bromeliads
catlist <- 
  bwgdet %>% 
  dplyr::select(bromeliad_id, cat, min, max) %>% 
  dplyr::group_by(bromeliad_id) %>% 
  dplyr::tally()
## Those with one category
bwgdet %>% 
  dplyr::filter(bromeliad_id %in% (catlist %>% 
                                     dplyr::filter(n == 1) %>% 
                                     dplyr::select(bromeliad_id) %>% 
                                     dplyr::pull())) %>% 
  dplyr::select(cat) %>% 
  dplyr::group_by(cat) %>% 
  tally()
### Those with two categories
bwgdet %>% 
  dplyr::filter(bromeliad_id %in% (catlist %>% 
                                     dplyr::filter(n == 2) %>% 
                                     dplyr::select(bromeliad_id) %>% 
                                     dplyr::pull()))%>% 
  dplyr::select(cat) %>% 
  dplyr::group_by(cat) %>% 
  tally()
### Those with three categories
bwgdet %>% 
  dplyr::filter(bromeliad_id %in% (catlist %>% 
                                     dplyr::filter(n == 3) %>% 
                                     dplyr::select(bromeliad_id) %>% 
                                     dplyr::pull())) %>% 
  dplyr::select(cat) %>% 
  dplyr::group_by(cat) %>% 
  tally()

# Get those bromeliads with just coarse and fine
goodones <- 
  bwgdet %>% 
  ## Bromeliads with two size categories
  dplyr::filter(bromeliad_id %in% (bwgdet %>% 
                                     dplyr::filter(bromeliad_id %in% (catlist %>% 
                                                                        dplyr::filter(n == 2) %>% 
                                                                        dplyr::select(bromeliad_id) %>% 
                                                                        dplyr::pull())) %>% 
                                     dplyr::select(bromeliad_id, cat) %>%
                                     dplyr::filter(cat != "all") %>% 
                                     dplyr::group_by(bromeliad_id) %>% 
                                     dplyr::tally() %>% 
                                     dplyr::filter(n == 2) %>% 
                                     dplyr::select(bromeliad_id) %>% 
                                     dplyr::pull())) 
## Can't find the code to estimate below 100, so pick overall value in fwdata
## do the difference to get FPOM with estimation
goodonesfw <- 
  fwdata$bromeliads %>% 
  #dplyr::filter(visit_id == 486) %>% 
  dplyr::select(visit_id, bromeliad_id, min_det, max_det, total_detritus) %>%
  dplyr::left_join(goodones %>% 
                     dplyr::select(-min, -max) %>% 
                     tidyr::pivot_wider(values_from = mass_g,
                                        names_from = cat),
                   by = c("bromeliad_id", "visit_id")) %>% 
  dplyr::mutate(fine_sub = total_detritus - coarse)

  


## Look for potential categories
## Extract bromeliads that have uneven grouping (i.e. something else than
## just fine and coarse)


### APPLY CORRECTION FOR THOSE NOT GOING TO 0



# Get ratio of coarse to fine organic matter in brastri bromeliads --------
# In brasil
brastridet <- 
  readr::read_csv(here::here("data",
                             "bromeliad_data.csv")) %>%
  ## Add detritus categories
  dplyr::mutate(coarse_det_g = loose_det_g + coarse_det_g) %>% 
  ## Select relevant columns
  dplyr::select(bromeliad_id, country, bromspecies,
                coarse_det_g, fine_det_g) %>% 
  ## Get detritus ratio
  dplyr::mutate(det_ratio = coarse_det_g/fine_det_g) %>% 
  dplyr::filter(!is.na(det_ratio)) 

mean(brastridet$det_ratio)
median(brastridet$det_ratio)

plotly::ggplotly(
  brastridet %>% 
    ggplot(aes(x = det_ratio)) + 
    geom_density(fill = "darkgoldenrod1"))

# In the BWG bromeliads




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
# about to months to develop to adulthood, no precise measurement
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
676.5

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
270 + 420

# Lampyridae
# waiting reply from lead author
## Vaz, S., L. F. L. Da Silveira, and S. P. Rosa. 2020. Morphology and 
## life cycle of a new species of Psilocladus Blanchard, 1846 
## (Coleoptera, Lampyridae, Psilocladinae), the first known 
## bromeliad-inhabiting firefly. Pap. Avulsos Zool. 60: e202060(s.i.).24.
# so far last instar to emergence is 5 months

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
                             60, 676.5, 730, 450,
                             270, 200, 144, 333,
                             150, 
                             (60 + 676.5 + 730 + 450 + 270 + 200 + 333) / 7,
                             (60.5 + 150) / 2))


# Get CPI
traits <- 
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
  dplyr::rowwise() %>% 
  dplyr::mutate(CPI_crisp = ifelse(is.na(CP1)|is.na(CP2)|is.na(CP3),
                                   NA, mean(c(rep(10, CP1), rep(40, CP2), 
                                              rep(match_CP3(ord = ord, subord = subord, 
                                                            family = family, genus = genus,
                                                            cp3_search = cp3_search), 
                                                  CP3)))))






# stoich data -------------------------------------------------------------
stoich <- 
  readr::read_csv(here::here("data",
                             "stoichiometry.csv")) %>% 
  dplyr::select(Site, Class, Order, Family, Genus,Sp) %>% 
  dplyr::distinct()
readr::write_csv(stoich,
                 here::here("angélica_stoich_species.csv"))




