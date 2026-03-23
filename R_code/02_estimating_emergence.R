# Estimating emergence of bromeliad insects

# Libraries
library(tidyverse)
library(here)


# Load and prep instar data -----------------------------------------------
# Data from Olivier's paper, some code is from there too!
## Dézerald, O., C. Leroy, B. Corbara, A. Dejean, S. Talaga, and R. Céréghino. 
## 2017. Environmental drivers of invertebrate population dynamics in Neotropical 
## tank bromeliads. Freshw Biol 62: 229–242.

## First let's make a named list with every organism of interest
species <- 
  c("Anopheles", "Bezzia", "Cerato",
    "Chiro", "Corethrella", "Culex",
    "Cyphon", "Limonidae", "WyeomyiiaSP1",
    "WyeomyiiaSP2")


## Read the .txt data for each of these sepcies
instars <- 
  species %>% 
  purrr::map(
    ~readr::read_tsv(here::here("data",
                                "oliviemergence",
                                paste0(.x,
                                       ".txt")))) %>% 
  ### rename columns for ease of manipulation
  purrr::map(
    ~.x %>% 
      dplyr::rename(graduation = dplyr::ends_with("grad"),
                    size_mm = dplyr::ends_with("mm"))) %>% 
  ### Convert from French to English decimal system
  purrr::map(
    ~.x %>% 
      dplyr::mutate(size_mm = stringr::str_replace(size_mm,
                                              ",",
                                              "."),
                    graduation = stringr::str_replace(size_mm,
                                                   ",",
                                                   "."))) %>% 
  ### Make columns numeric
  purrr::map(
    ~.x %>% 
      dplyr::mutate(size_mm = as.numeric(size_mm)))

## Set names because maybe it is useful
names(instars) <- 
  species

## Reuse the code from Olivier to establish instars
## Gonna have to do this by hand
## Anopheles
instars$Anopheles <- 
  instars$Anopheles %>% 
  dplyr::mutate(instar = ifelse(size_mm < 0.18, 1,
                                ifelse(size_mm >= 0.18 & size_mm <= 0.26, 2,
                                       ifelse(size_mm >= 0.36 & size_mm <= 0.44, 3,
                                              ifelse(size_mm > 0.44, 4, NA)))))
## Bezzia
instars$Bezzia <- 
  instars$Bezzia %>% 
  dplyr::mutate(instar = ifelse(size_mm <=0.16, 1,
                                ifelse(size_mm > 0.16 & size_mm < 0.22, 2,
                                       ifelse(size_mm >= 0.22 & size_mm <=0.30, 3,
                                              ifelse(size_mm > 0.30, 4, NA)))))

## Cerato
instars$Cerato <- 
  instars$Cerato %>% 
  dplyr::mutate(instar = ifelse(size_mm < 0.20, 1,
                                ifelse(size_mm >= 0.20 & size_mm <= 0.28, 2,
                                       ifelse(size_mm > 0.28, 3, NA))))

## Chiro
instars$Chiro <- 
  instars$Chiro %>% 
  dplyr::mutate(instar = ifelse(size_mm < 0.18, 1,
                                ifelse(size_mm >= 0.18 & size_mm <= 0.26, 2,
                                       ifelse(size_mm > 0.26, 3, NA))))

## Corethrella
instars$Corethrella <- 
  instars$Corethrella %>% 
  dplyr::mutate(instar = ifelse(size_mm <= 0.22, 1,
                                ifelse(size_mm > 0.22 & size_mm <= 0.38, 2,
                                       ifelse(size_mm > 0.38 & size_mm <= 0.58, 3,
                                              ifelse(size_mm > 0.58, 4, NA)))))
  
 
## Culex
instars$Culex <- 
  instars$Culex %>% 
  dplyr::mutate(instar = ifelse(size_mm < 0.38, 1,
                                ifelse(size_mm >= 0.38 & size_mm <= 0.54, 2,
                                       ifelse(size_mm > 0.54 & size_mm <= 0.66, 3,
                                              ifelse(size_mm >= 0.66, 4, NA)))))
 
  
## Cyphon  
instars$Cyphon <- 
  instars$Cyphon %>% 
  dplyr::mutate(instar = ifelse(size_mm <= 0.18, 1,
                                ifelse(size_mm > 0.18 & size_mm <= 0.28, 2,
                                       ifelse(size_mm > 0.28 & size_mm <= 0.38, 3,
                                              ifelse(size_mm > 0.38 & size_mm <= 0.48, 4,
                                                     ifelse(size_mm > 0.48 & size_mm <= 0.58, 5,
                                                            ifelse(size_mm > 0.58 & size_mm <= 0.68, 6,
                                                                   ifelse(size_mm > 0.68 & size_mm  <= 0.78, 7,
                                                                          ifelse(size_mm > 0.78 & size_mm <= 0.88,8,
                                                                                 ifelse(size_mm > 0.88 & size_mm <= 0.98, 9,
                                                                                        ifelse(size_mm > 0.98, 10, NA)))))))))))

## Limoniidae
instars$Limonidae <- 
  instars$Limonidae %>% 
  dplyr::mutate(instar = ifelse(size_mm <= 0.14, 1,
                                ifelse(size_mm > 0.14 & size_mm <=0.22, 2,
                                       ifelse(size_mm > 0.22 & size_mm <= 0.48, 3,
                                              ifelse(size_mm > 0.48, 4, NA)))))

## Wyeomyiia species 1
instars$WyeomyiiaSP1 <- 
  instars$WyeomyiiaSP1 %>% 
  dplyr::mutate(instar = ifelse(size_mm < 0.36, 1,
                                ifelse(size_mm >= 0.36 & size_mm < 0.50, 2,
                                       ifelse(size_mm >= 0.50 & size_mm <= 0.68, 3,
                                              ifelse(size_mm > 0.68, 4, NA)))))

## Wyeomyiia species 2
instars$WyeomyiiaSP2 <- 
  instars$WyeomyiiaSP2 %>% 
  dplyr::mutate(instar = ifelse(size_mm < 0.30, 1,
                                ifelse(size_mm >= 0.30 & size_mm <= 0.42, 2,
                                       ifelse(size_mm > 0.42 & size_mm <= 0.58, 3,
                                              ifelse(size_mm > 0.58, 4, NA)))))
  

## Now flatten list and keep the column that we want only
instars <- 
  instars %>% 
  purrr::list_rbind(names_to = "species") %>% 
  dplyr::select(date = Date, nplant = NPlantID, species,
                instar, graduation, size_mm)

## Make sure no instar is NA
stopifnot(!is.na(instars$instar))



# Check pupa data -----------------------------------------------------------
# Read data
pupae <- 
  readr::read_csv(here::here("data",
                             "oliviemergence",
                             "pupae.csv"))

## Check columns
colnames(pupae)

# Combine columns
pupae <- 
  pupae %>% 
  ## Combine culicid unknwon columns
  ## do it row by row
  dplyr::rowwise() %>% 
  dplyr::mutate(culicid_a = sum(dplyr::c_across(dplyr::contains("Culicidae A")),
                                na.rm = TRUE),
                culicid_b = sum(dplyr::c_across(dplyr::contains("Culicidae B")),
                                na.rm = TRUE),
                petite_jaune = sum(dplyr::c_across(dplyr::contains("bifide")),
                             na.rm = TRUE),
                petite_blanche = sum(dplyr::c_across(dplyr::contains("syphon")),
                              na.rm = TRUE)) %>% 
  ## remove row grouping
  dplyr::ungroup() %>% 
  ## Remove old ones
  dplyr::select(-Date,
                -"Plant ID",
                -dplyr::contains("Culicidae A"),
                -dplyr::contains("Culicidae B"),
                -dplyr::contains("bifide"),
                -dplyr::contains("syphon"),
                -"Cyphon adulte") %>% 
  ## rename cyphon and bromeliad columns
  dplyr::rename(cyphon = "Cyphon nymphe",
                nplant = NplantID) %>% 
  ## remove NAs
  dplyr::filter(!is.na(nplant)) %>% 
  ## Add who was sampled when
  dplyr::left_join(instars %>% 
                     dplyr::select(nplant, date) %>% 
                     dplyr::distinct(),
                   by = "nplant")


## Find which species is which pupa comparing to Figure 4 of Dézerald et al. 
## 2017 10.1111/fwb.12862
## Using Dynamics - Final Graphes2.R from Olivier, looks like wyeo1 is aphobema
## and wyeo2 is lamellata
pupae %>%
  ## Sum by all columns by sampling date
  dplyr::group_by(date, nplant) %>%
  dplyr::summarise(dplyr::across(c(Corethrella:culicid_b),
                                 ~sum(.))) %>% 
  dplyr::ungroup() %>% 
  ## Convert date to date format to have in the right order
  dplyr::mutate(date = lubridate::dmy(date)) %>%
  ## order by ascending date
  dplyr::arrange(date)
  
  
## My five unknown nymphs
## wyeomyiia lamellata/sp2 (nothing until second sample only) cerato (3 samples w pupae)
## Culicidae C is Anopheles (nothing until 7th sample)
## Culicidae A is Wywomyia aphobema/sp1, first two samples have pupae
## Culicidae B is culex, nothing until 4th sample
## petite jaune 
## petite blanche

## Finish cleaning the data
pupae <- 
  pupae %>% 
  ## Remove petites until because we will never know what they are
  dplyr::select(-petite_jaune,
                -petite_blanche) %>%
  ### Rename columns to match those of instars
  dplyr::rename(Anopheles = "Culicidae C",
                WyeomyiiaSP1 = culicid_a,
                Culex = culicid_b,
                Cyphon = cyphon,
                Limonidae = Limo) %>% 
  ### Pivot longer to match instars format
  tidyr::pivot_longer(cols = Corethrella:Culex,
                      names_to = "species",
                      values_to = "n") %>%
  ### Remove rows with zeroes
  dplyr::filter(n > 0) %>%
  ### Expand rows by number of pupae
  tidyr::uncount(n) %>% 
  ## Add instar column
  dplyr::left_join(instars %>% 
                     ## Select columns of interest
                     dplyr::select(species, instar) %>% 
                     ## get unique instars per species
                     dplyr::distinct() %>%
                     ## only keep the highest one
                     dplyr::group_by(species) %>% 
                     dplyr::filter(instar == max(instar)) %>% 
                     dplyr::ungroup() %>% 
                     ## add one to match pupal instar
                     dplyr::mutate(instar = instar +1)                     ,
                   by = c("species"))

# Combine instar and pupa data 
## Compared to the overall number of specimens collected, there are barely any pupae (~1%)
## So lets do full curves
nrow(pupae)
nrow(instars)

# Indiviudal level
instapup <- 
  instars %>% 
  ## Add pupae
  dplyr::bind_rows(pupae) %>% 
  ## Fix these typos
  dplyr::mutate(species = ifelse(species == "Limonidae",
                                 "Limoniidae", species)) %>% 
  dplyr::mutate(species = ifelse(species == "WyeomyiiaSP1",
                                 "WyeomyiaSP1", species)) %>% 
  ## Remove cerato and wyeomyiia sp2
  dplyr::filter(!(species %in% c("Cerato", 
                                 "WyeomyiiaSP2"))) 
  

# Pooled per stage
instapup_pooled <- 
  instapup %>% 
  dplyr::select(-date, -nplant, - graduation, -size_mm) %>% 
  dplyr::group_by(species, instar) %>%
  dplyr::tally()

# Pooled per stage, pupae with last stage
instapup_pooled_combined <- 
  instapup %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(instar = ifelse(instar ==max(instar),
                                max(instar) - 1, instar)) %>% 
  dplyr::select(-date, -nplant, - graduation, -size_mm) %>% 
  dplyr::group_by(species, instar) %>%
  dplyr::tally()

# Pool by cohort
## Will be hard to do, data not gettable for most species

# Pupae are so rare because it is probably a very short stage

## So let's go without pupae, keep this code here in case of


# Combine data per instar without pupae --------------------------------------
# Indiviudal level
instars_pooled <- 
  instars %>% 
  ## Fix these typos
  dplyr::mutate(species = ifelse(species == "Limonidae",
                                 "Limoniidae", species)) %>% 
  dplyr::mutate(species = ifelse(species == "WyeomyiiaSP1",
                                 "WyeomyiaSP1", species))  %>% 
  ## Pool per instar
  dplyr::select(-date, -nplant, - graduation, -size_mm) %>% 
  dplyr::group_by(species, instar) %>%
  dplyr::tally()
  

# Assess survival rate at each level --------------------------------------
# Get number of specimens across each level
instars_pooled %>% 
  ggplot(aes(x = instar,
             y = n)) +
  geom_smooth() +
  facet_wrap(~species, 
             scale = "free_x")

# Remove those that have a non-decreasing number of specimes
instars_pooled <- 
  instars_pooled %>% 
  ## Remove cerato and wyeomyiia sp2
  dplyr::filter(!(species %in% c("Anopheles", 
                                 "Corethrella",
                                 "Culex",
                                 "WyeomyiiaSP2"))) 

# Estimate relationships to one more instar, ie adulthood -----------------
## Get biomass for the different instars
larvalmass <- 
  ## Read data
  readr::read_table(here::here("data",
                               "oliviemergence",
                               "weight.txt"),
                    ## because of the commas (French system)
                    ## read all columns as character
                    col_types = cols(.default = col_character())) %>% 
  ## convert commas to dots and convert to numeric
  dplyr::mutate(dplyr::across(dplyr::contains("Weight"),
                              ~stringr::str_replace(.,
                                                    ",",
                                                    "."))) %>%
  dplyr::mutate(dplyr::across(dplyr::contains(c("Weight", "Nb")),
                              ~as.numeric(.))) %>%
  ## Now try to understand if larval weight is just one larva or all of a sample
  dplyr::mutate(test = AluLarvaeWeight - AluWeight,
                test2 = (AluLarvaeWeight - AluWeight)/Nblarvae) %>% 
  dplyr::mutate(test = test - LarvaeWeight.mg.,
                test2 = test2 - LarvaeWeight.mg.) %>% 
  ## Ok so it looks like it is the average mass of larvae
  ## so select and rename
  dplyr::select(species = Larval_Stage,
                avg_mass_mg = LarvaeWeight.mg.) %>% 
  ## Now we just need to split by instar, so get last character in the string
  dplyr::mutate(instar = as.numeric(stringr::str_sub(species,
                                          -1))) %>% 
  ## Convert those with instar 10 (so 0 with command above) to 10
  dplyr::mutate(instar = ifelse(instar == 0,
                                10, instar)) %>% 
  ## Now convert species to match our naming system
  dplyr::mutate(species = ifelse(stringr::str_detect(string = species,
                                                     pattern = "Culex"),
                                 "Culex", ifelse(stringr::str_detect(string = species,
                                                                     pattern = "Cerato"),
                                                 "Cerato", ifelse(stringr::str_detect(string = species,
                                                                                         pattern = "Bez"),
                                                                     "Bezzia", ifelse(stringr::str_detect(string = species,
                                                                                                          pattern = "Chiro"),
                                                                                      "Chiro", ifelse(stringr::str_detect(string = species,
                                                                                                                          pattern = "Cyp"),
                                                                                                      "Cyphon", ifelse(stringr::str_detect(string = species,
                                                                                                                                            pattern = "Core"),
                                                                                                                        "Corethrella", ifelse(stringr::str_detect(string = species,
                                                                                                                                                                  pattern = "Limo"),
                                                                                                                                              "Limoniidae", ifelse(stringr::str_detect(string = species,
                                                                                                                                                                                        pattern = "SP1"),
                                                                                                                                                                    "WyeomyiaSP1", species))))))))) %>% 
  ## Finally keep only the species we want
  dplyr::filter(species %in% unique(instars_pooled$species))

# Stick to the pooled data and estimate biomass
instars_pooled <- 
  instars_pooled %>% 
  dplyr::left_join(larvalmass,
                   by = c("species", "instar")) %>%
  ## Get estimated biomass
  dplyr::mutate(est_biomass_mg= n * avg_mass_mg)

# Check the plot
instars_pooled %>% 
  ggplot(aes(x = instar,
             y = est_biomass_mg)) +
  geom_smooth() +
  facet_wrap(~species, 
             scale = "free_x")

# Check log log
instars_pooled %>% 
  ggplot(aes(x = log(instar),
             y = log(est_biomass_mg))) +
  geom_smooth() +
  facet_wrap(~species, 
             scale = "free_x")


# Get survivorship of whole cohort to last instar ------------------------------------
## Get average survival rate for all species
s_rates <- 
  instars_pooled %>% 
  ## Split per species
  dplyr::group_by(species) %>% 
  dplyr::group_split() %>% 
  ## Divide each value wiht the one right after in the vector (i.e. one instar above)
  purrr::map(~ .x %>% 
               ## Make sure everything is well sorted
               dplyr::arrange(species, instar) %>% 
               dplyr::group_by(species) %>% 
               ## Get survivorship rate
               dplyr::mutate(survivorship_rate = dplyr::lead(n) / n) %>% 
               ## Make valeues above 1 NAs
               dplyr::mutate(survivorship_rate = ifelse(survivorship_rate >= 1,
                                                        NA, survivorship_rate)) %>% 
               ## Remove NAs
               dplyr::filter(!is.na(survivorship_rate))) %>% 
  ## Join back together
  purrr::list_rbind() %>% 
  ## Get average survival rate between two instars
  dplyr::select(species, survivorship_rate) %>% 
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~mean(.)))

# Get how many would get to last instar + 1 (adults)
ol_probem <- 
  instars_pooled %>%
  ## Select species and individuals
  dplyr::select(species, n, est_biomass_mg) %>% 
  ## Get how many individuals overall in the species
  dplyr::summarise(dplyr::across(c(n, 
                                   est_biomass_mg),
                                 ~sum(.))) %>% 
  ## Add expected number of adults
  dplyr::left_join(instars_pooled %>% 
                     ### Group by species
                     dplyr::group_by(species) %>% 
                     ### Select relevant columns
                     dplyr::select(species, instar, n, avg_mass_mg) %>% 
                     ### Filter to last instar
                     dplyr::filter(instar == max(instar)) %>% 
                     ### Join survival rates
                     dplyr::left_join(s_rates,
                                      by = "species") %>% 
                     ### Get how many individual would make it to one more instar
                     ### Using the estimated rates from above
                     dplyr::mutate(n = n * survivorship_rate) %>% 
                     ### Get biomass lost through these adults
                     dplyr::mutate(biomass_lost_mg = n * avg_mass_mg) %>% 
                     ### Select relevant columns
                     dplyr::select(species, n_adults = n, biomass_lost_mg),
                   by = "species") %>% 
  ## Get probability of becoming an adult
  dplyr::mutate(probem = n_adults / n,
                prop_biomass_lost = biomass_lost_mg / est_biomass_mg) %>%
  ## Select relevant columns
  dplyr::select(species, probem, prop_biomass_lost)



# Add to BWG data, up to order level, then see whats up -------------------
## Add to BWG data
bwg_probem <- 
  readr::read_rds(here::here("data",
                             "scalewebs",
                             "fw_database.rds")) %>% 
  ## Keep only the trait element that has all the species
  purrr::pluck("traits") %>% 
  ## Keep relevant columns
  dplyr::select(bwg_name, class, ord,
                family, genus) %>% 
  ## Add the column, first saying that non-insects will not emerge
  dplyr::mutate(probem = ifelse(class != "Insecta",
                                0, NA),
                prop_biomass_lost = ifelse(class != "Insecta",
                                         0, NA)) %>% 
  ## Some classes are NA for non-insects
  dplyr::mutate(probem = ifelse(is.na(class),
                                0, probem),
                prop_biomass_lost = ifelse(is.na(class),
                                         0, prop_biomass_lost)) %>%
  ## Add that Dytiscidae, Elmidae, Hydrophilidae, Veliidae (kinda) are entirely aquatic
  dplyr::mutate(probem = ifelse(family %in% c("Dytiscidae", "Elmidae",
                                              "Hydrophilidae", "Veliidae", "Vellidae",
                                              "Gerridae", "Mesoveliidae"),
                                0, probem),
                prop_biomass_lost = ifelse(family %in% c("Dytiscidae", "Elmidae",
                                                       "Hydrophilidae", "Veliidae", "Vellidae",
                                                       "Gerridae", "Mesoveliidae"),
                                         0, prop_biomass_lost)) %>% 
  ## Add 0 to the semi-aquatic/terrestrial families, hard to say
  dplyr::mutate(probem = ifelse(!is.na(family) & family %in% 
                                  c("Carabidae", "Curculionidae", 
                                    "Elateridae", "Enicocephalidae",
                                    "Tineidae", "Scutelleridae"),
                                0, probem),
                prop_biomass_lost = ifelse(!is.na(family) & family %in% 
                                             c("Carabidae", "Curculionidae", 
                                               "Elateridae", "Enicocephalidae",
                                               "Tineidae", "Scutelleridae"),
                                           0, prop_biomass_lost))

# Start at genus level
bwg_probem <- 
  bwg_probem %>% 
  ## Add genus level info
  dplyr::left_join(ol_probem %>%
                     ### Genus
                     dplyr::rename(genus = species,
                                   probem_genus = probem,
                                   prop_biomass_lost_genus = prop_biomass_lost) %>%
                     ### Fix Wyeomyia
                     dplyr::mutate(genus = ifelse(genus == "WyeomyiaSP1",
                                                  "Wyeomyia", genus)),
                   by = "genus") %>%
  ## Fill the NAs
  dplyr::mutate(probem = ifelse(is.na(probem),
                                probem_genus, probem),
                prop_biomass_lost = ifelse(is.na(prop_biomass_lost),
                                         prop_biomass_lost_genus, prop_biomass_lost)) %>%
  ## Remove column
  dplyr::select(-probem_genus, -prop_biomass_lost_genus)

# Do family level
bwg_probem <- 
  bwg_probem %>% 
  ## Add family level info
  dplyr::left_join(ol_probem %>%
                     ### Family
                     dplyr::mutate(family = species,
                                   probem_family = probem,
                                   prop_biomass_lost_family = prop_biomass_lost) %>%
                     ### Adapt taxonomy to new level
                     dplyr::mutate(family =  ifelse(family == "WyeomyiaSP1",
                                                    "Culicidae", ifelse(family %in% c("Bezzia", "Cerato"),
                                                                        "Ceratopogonidae", ifelse(family == "Chiro",
                                                                                                  "Chironomidae", ifelse(family == "Cyphon",
                                                                                                                         "Scirtidae", family))))) %>%
                     ### Add line for Tipulidae, which is essentially Limoniidae
                     dplyr::bind_rows(tibble::tibble(family = "Tipulidae",
                                                     probem_family = ol_probem %>%
                                                       dplyr::filter(species == "Limoniidae") %>%
                                                       dplyr::pull(probem),
                                                     prop_biomass_lost_family = ol_probem %>%
                                                       dplyr::filter(species == "Limoniidae") %>%
                                                       dplyr::pull(prop_biomass_lost))) %>%
                     ### Summarise per family
                     dplyr::group_by(family) %>% 
                     dplyr::summarise(dplyr::across(c(probem_family, prop_biomass_lost_family),
                                                    ~mean(.))),
                   by = "family") %>%
  
  ## Fill the NAs
  dplyr::mutate(probem = ifelse(is.na(probem),
                                probem_family, probem),
                prop_biomass_lost = ifelse(is.na(prop_biomass_lost),
                                           prop_biomass_lost_family, prop_biomass_lost)) %>%
  ## Rename column
  dplyr::select(-probem_family, -prop_biomass_lost_family)


# Lit review for those other species without estimate --------------------------------------
# Pull out known families that are still NA
bwg_probem %>% 
  dplyr::filter(is.na(probem)) %>% 
  dplyr::select(family) %>% 
  dplyr::distinct()

# Odonata
## Lawton, J.H. (1970) A population study on larvae of the damselfly Pyrrhosoma 
## nymphula (Sulzer) (Odonata: Zygoptera). Hydrobiologia 36, 33–52. 
## https://doi.org/10.1007/BF00751281
## Combined the values obtained in Lawton Fig. 3 above with average biomass for each instars
## found here in Fig. 4
## Waringer, J.A. and Humpesch, U.H. (1984), Embryonic development, larval growth 
## and life cycle of Coenagrion puella (Odonata: Zygoptera) from an Austrian pond. 
## Freshwater Biology, 14: 385-399
## https://doi.org/10.1111%2Fj.1365-2427.1984.tb00162.x

## Make data
odonates <- 
  tibble::tibble(instar = c(1:11),
                 n = c(550, 175, 355, 165, 195, 87, 60, 44, 41, 35, 29),
                 avg_body_mass_mg = c(0.02, 0.03, 0.055, 0.09, 0.16, 0.26, 0.5, 0.9, 1.7, 3.2, 6)) %>% 
  ### Get biomass at each trophic level
  dplyr::mutate(biomass_mg = n * avg_body_mass_mg) %>% 
  ### Sum everything
  dplyr::summarise(dplyr::across(c(n, biomass_mg),
                                 ~sum(.))) %>% 
  ## Now get proportion based on values from the two papers
  dplyr::mutate(n = 23/n, biomass_mg = 138/biomass_mg)

## Add to data
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(!is.na(ord) & ord == "Odonata",
                                odonates$n, probem),
                prop_biomass_lost = ifelse(!is.na(ord) & ord == "Odonata",
                                           odonates$biomass_mg, prop_biomass_lost))




# Calamoceratidae
## No data on family
## For Trichoptera as a whole:
## Otto, C. (1974). Growth and Energetics in a Larval Population of Potamophylax 
## Cingulatus (Steph.) (Trichoptera) in a South Swedish Stream. 
## Journal of Animal Ecology, vol. 43, no. 2, pp. 339–61.
## https://doi-org.proxy3.library.mcgill.ca/10.2307/3369. 
## Life table partial and in Sweden but nothing else

## Using figure 2
## 44/785 individuals were adults
## 825/3404.42 biomass measured was in adults
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(!is.na(ord) & ord == "Trichoptera",
                                44/785, probem),
                prop_biomass_lost = ifelse(!is.na(ord) & ord == "Trichoptera",
                                           825/3404.42, prop_biomass_lost))


# Lampyridae
## Lawrent L. Buschman, Larval Development and Its Photoperiodic Control in 
## the Firefly Pyractomena lucifera (Coleoptera: Lampyridae), Annals of the 
## Entomological Society of America, Volume 81, Issue 1, 1 January 1988, 
## Pages 82–90, https://doi.org/10.1093/aesa/81.1.82
## Using mass for each instar obtained in rear studies (fig 1A),
## combined with proportion of individuals in each instar collected from the field (Fig 3)
lampyridae <- 
  tibble::tibble(instars = c(1:5),
                 avg_body_mass_mg = c(0.65, 1.75, 5.4, 14, 47.5),
                 n = c(17, 17, 14, 17, 24)) %>% 
  ### Get biomass at each trophic level
  dplyr::mutate(biomass_mg = n * avg_body_mass_mg) %>% 
  ### Sum everything
  dplyr::summarise(dplyr::across(c(n, biomass_mg),
                                 ~sum(.))) %>% 
  ## Now get proportion based on values from the paper
  dplyr::mutate(n = 0.1, biomass_mg = 465/biomass_mg)

## Add to data
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(!is.na(family) & family == "Lampyridae",
                                lampyridae$n, probem),
                prop_biomass_lost = ifelse(!is.na(family) & family == "Lampyridae",
                                           lampyridae$biomass_mg, prop_biomass_lost))




# Ptylodactylidae
## LeSage, Laurent, & Harper, P P. (1976). Notes on the Life History of the 
## Toed-Winged Beetle Anchytarsus bicolor (Melsheimer) (Coleoptera: Ptilodactylidae). 
## The Coleopterists’ Bulletin, 30(3), 233--238. https://doi.org/10.5962/p.371986
## But cannot be combined with, because pronotal length vs. full body length 
##Benke, A. C., A. D. Huryn, L. A. Smock, and J. B. Wallace. 1999. Length-mass 
## relationships for freshwater macroinvertebrates in North America with particular 
## reference to the southeastern United States. Journal of the North American 
## Benthological Society 18: 308–343.

## Apply average of other coleopterans
### Get average for coleopterans
avg_coleo <- 
  bwg_probem %>% 
  ### filter Coleopterans
  dplyr::filter(ord == "Coleoptera") %>% 
  ### Remove NAs
  dplyr::filter(!is.na(probem))%>%
  ## Keep columns of interest
  dplyr::select(probem, prop_biomass_lost) %>% 
  ## Get unique values
  dplyr::distinct() %>% 
  ## Get mean
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~mean(.)))
## Add it to data
bwg_probem <- 
  bwg_probem %>% 
  ### Apply values
  dplyr::mutate(probem = ifelse(!is.na(family) & family == "Ptilodactylidae",
                                avg_coleo$probem, probem),
                prop_biomass_lost = ifelse(!is.na(family) & family == "Ptilodactylidae",
                                           avg_coleo$prop_biomass_lost, prop_biomass_lost))

## Add the same value to all unidentified Coleoptera
bwg_probem <- 
  bwg_probem %>% 
  ### Apply values
  dplyr::mutate(probem = ifelse(ord == "Coleoptera" & is.na(probem),
                                avg_coleo$probem, probem),
                prop_biomass_lost = ifelse(ord == "Coleoptera" & is.na(prop_biomass_lost),
                                           avg_coleo$prop_biomass_lost, prop_biomass_lost))

# Other Dipterans without references
## NB Ephydridae has one study but in lab and survical very high (90%+)
## Ugine TA, Sanderson JP, Wraight SP. Developmental times and life tables for 
## shore flies, Scatella tenuicosta (Diptera: Ephydridae), at three temperatures. 
## Environ Entomol. 2007 Oct;36(5):989-97. doi: 10.1603/0046-225x(2007)36[989:dtaltf]2.0.
## Stratomyidae only controlled studies on Hermetia illucens
## Syrphidae only studies in controlled lab environments
## amdem, M. M., & Otomo, P. V. (2023). Developmental performance of Eristalis 
## tenax larvae (Diptera: Syrphidae): Influence of growth media and yeast addition 
## during captive rearing. Journal of Experimental Zoology Part A: Ecological and 
## Integrative Physiology, 339, 503–513. https://doi.org/10.1002/jez.2696
## aff. Drosophilidae because ID uncertain
## Anisopodidae has https://doi.org/10.1016/S0043-1354(96)00236-9 for distribution
## of body size cohort, but no allometric equation to apply afterwars
## Apply average of other coleopterans
### Get average for coleopterans
avg_dipte <- 
  bwg_probem %>% 
  ### Filter dipterans
  dplyr::filter(ord == "Diptera") %>% 
  ### Remove NAs
  dplyr::filter(!is.na(probem))%>%
  ## Keep columns of interest
  dplyr::select(probem, prop_biomass_lost) %>% 
  ## Get unique values
  dplyr::distinct() %>% 
  ## Get mean
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~mean(.)))
## Add it to data
bwg_probem <- 
  bwg_probem %>% 
  ### Apply values
  dplyr::mutate(probem = ifelse(!is.na(family) & family %in% c("Axymyiidae", "Canacidae",
                                                               "Cecidomyiidae", "Sciaridae",
                                                               "Corethrellidae", "Dolichopodidae",
                                                               "Empididae", "Ephydridae", "Muscidae",
                                                               "Periscelididae", "Phoridae",
                                                               "Psychodidae", "Sarchophagidae", ## typo
                                                               "Scatopsidae", "Sphaeroceridae",
                                                               "Stratiomyidae", "Tabanidae",
                                                               "Thaumaleidae", "Aulacigastridae",
                                                               "Dixidae", "Bibionidae", "Anisopodidae",
                                                               "Therevidae", "aff. Drosophilidae",
                                                               "Syrphidae"),
                                avg_dipte$probem, probem),
                prop_biomass_lost = ifelse(!is.na(family) & family %in% c("Axymyiidae", "Canacidae",
                                                                          "Cecidomyiidae", "Sciaridae",
                                                                          "Corethrellidae", "Dolichopodidae",
                                                                          "Empididae", "Ephydridae", "Muscidae",
                                                                          "Periscelididae", "Phoridae",
                                                                          "Psychodidae", "Sarchophagidae", ## typo
                                                                          "Scatopsidae", "Sphaeroceridae",
                                                                          "Stratiomyidae", "Tabanidae",
                                                                          "Thaumaleidae", "Aulacigastridae",
                                                                          "Dixidae", "Bibionidae", "Anisopodidae",
                                                                          "Therevidae", "aff. Drosophilidae",
                                                                          "Syrphidae"),
                                          avg_dipte$prop_biomass_lost, prop_biomass_lost))   

## Add the same value to all unidentified Diptera
bwg_probem <- 
  bwg_probem %>% 
  ### Apply values
  dplyr::mutate(probem = ifelse(ord == "Diptera" & is.na(probem),
                                avg_dipte$probem, probem),
                prop_biomass_lost = ifelse(ord == "Diptera" & is.na(prop_biomass_lost),
                                           avg_dipte$prop_biomass_lost, prop_biomass_lost))

# Now what is left is Hemiptera, Ephemeroptera, Lepidoptera, all unided, put 0
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(is.na(probem),
                                0, probem),
                prop_biomass_lost = ifelse(is.na(prop_biomass_lost),
                                           0, prop_biomass_lost))

# save for email
readr::write_csv(bwg_probem,
                 here::here("data",
                            "bwg_probem.csv"))




  