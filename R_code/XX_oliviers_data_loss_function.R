# Computing a loss function from Olivier's data

# Libraries
library(tidyverse)
library(here)


# Load and prep instar data -----------------------------------------------
## Reuse Olivier's code for this, we have the instars and the dates
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



# Add pupa data -----------------------------------------------------------
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


##  Find which species is which pupa comparing to Figure 4 of Dézerald et al. 
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


# Combine instar and pupa data --------------------------------------------
## Compared to initial number of specimen, there are barely any pupae (~1%)
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


# Assess survival rate at each level --------------------------------------
# Get number of specimens across each level
instapup_pooled %>% 
  ggplot(aes(x = instar,
             y = n)) +
  geom_smooth() +
  facet_wrap(~species, 
             scale = "free_x")

# Combine last two instars
instapup_pooled_combined %>% 
  ggplot(aes(x = instar,
             y = n)) +
  geom_smooth() +
  facet_wrap(~species, 
             scale = "free")

# So compile last instar and pupae together

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
                                                                     pattern = "Ano"),
                                                 "Anopheles", ifelse(stringr::str_detect(string = species,
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
  dplyr::filter(species %in% unique(instapup_pooled_combined$species))

# Stick to the pooled data and estimate biomass
instapup_pooled_combined <- 
  instapup_pooled_combined %>% 
  dplyr::left_join(larvalmass,
                   by = c("species", "instar")) %>%
  ## Get estimated biomass
  dplyr::mutate(est_biomass_mg= n * avg_mass_mg)

# Check the plot
instapup_pooled_combined %>% 
  ggplot(aes(x = instar,
             y = est_biomass_mg)) +
  geom_smooth() +
  facet_wrap(~species, 
             scale = "free_x")

# Check log log
instapup_pooled_combined %>% 
  ggplot(aes(x = log(instar),
             y = log(est_biomass_mg))) +
  geom_smooth() +
  facet_wrap(~species, 
             scale = "free_x")


# Get survivorship of whole cohort to last instar ------------------------------------
## Get average survival rate for all species
s_rates <- 
  instapup_pooled_combined %>% 
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
  instapup_pooled_combined %>%
  ## Select species and individuals
  dplyr::select(species,n) %>% 
  ## Get how many individuals overall in the species
  dplyr::summarise(dplyr::across(n,
                                 ~sum(.))) %>% 
  ## Add expected number of adults
  dplyr::left_join(instapup_pooled_combined %>% 
                     ###G roup by species
                     dplyr::group_by(species) %>% 
                     ### Select relevant columns
                     dplyr::select(species, instar, n) %>% 
                     ### Filter to last instar
                     dplyr::filter(instar == max(instar)) %>% 
                     ### Join survival rates
                     dplyr::left_join(s_rates,
                                      by = "species") %>% 
                     ### Get how many individual would make it to one more instar
                     ### Using the estimated rates from above
                     dplyr::mutate( n = n * survivorship_rate) %>% 
                     ### Select relevant columns
                     dplyr::select(species, n_adults = n),
                   by = "species") %>% 
  ## Get probability of becoming an adult
  dplyr::mutate(probem = n_adults / n) %>%
  ## Select relevant columns
  dplyr::select(species, probem) %>% 
  ## Remove NA
  dplyr::filter(!is.na(probem))



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
                                0, NA)) %>% 
  ## Some classes are NA for non-insects
  dplyr::mutate(probem = ifelse(is.na(class),
                                0, probem)) %>%
  ## Add that Dytiscidae, Elmidae, Hydrophilidae, Veliidae (kinda) are entirely aquatic
  dplyr::mutate(probem = ifelse(family %in% c("Dytiscidae", "Elmidae",
                                              "Hydrophilidae", "Veliidae", "Vellidae",
                                              "Gerridae", "Mesoveliidae"),
                                0, probem))

# Start at genus level
bwg_probem <- 
  bwg_probem %>% 
  ## Add genus level info
  dplyr::left_join(ol_probem %>%
                     ### Genus
                     dplyr::rename(genus = species,
                                   probem_genus = probem) %>%
                     ### Fix Wyeomyia
                     dplyr::mutate(genus = ifelse(genus == "WyeomyiaSP1",
                                                  "Wyeomyia", genus)),
                   by = "genus") %>%
  ## Fill the NAs
  dplyr::mutate(probem = ifelse(is.na(probem),
                                probem_genus, probem)) %>%
  ## Remove column
  dplyr::select(-probem_genus)

# Do family level
bwg_probem <- 
  bwg_probem %>% 
  ## Add family level info
  dplyr::left_join(ol_probem %>%
                     ### Family
                     dplyr::mutate(family = species,
                                   probem_family = probem) %>%
                     ### Adapt taxonomy to new level
                     dplyr::mutate(family =  ifelse(family %in% c("Anopheles", "Wyeomyia",
                                                                   "Culex"),
                                                    "Culicidae", ifelse(family == "Bezzia",
                                                                        "Ceratopogonidae", ifelse(family == "Chiro",
                                                                                                  "Chironomidae", ifelse(family == "Cyphon",
                                                                                                                         "Scirtidae", family))))) %>%
                     ### Add line for Tipulidae, which is essentially Limoniidae
                     dplyr::bind_rows(tibble::tibble(family = "Tipulidae",
                                                     probem_family = ol_probem %>%
                                                       dplyr::filter(species == "Limoniidae") %>%
                                                       dplyr::pull(probem))) %>%
                     ### Summarise per family
                     dplyr::group_by(family) %>% 
                     dplyr::summarise(dplyr::across(probem_family,
                                                    ~mean(.))),
                   by = "family") %>%
  
  ## Fill the NAs
  dplyr::mutate(probem = ifelse(is.na(probem),
                                probem_family, probem)) %>%
  ## Rename column
  dplyr::select(-probem_family)

# Do order level
bwg_probem <- 
  bwg_probem %>% 
  ## Add family level info
  dplyr::left_join(ol_probem %>%
                     ### Family
                     dplyr::mutate(ord = species,
                                   probem_ord = probem) %>%
                     ### Adapt taxonomy to new level
                     dplyr::mutate(ord =  ifelse(ord == "Cyphon",
                                                    "Coleoptera", "Diptera")) %>%
                     ### Summarise per order
                     dplyr::group_by(ord) %>% 
                     dplyr::summarise(dplyr::across(probem_ord,
                                                    ~mean(.))),
                   by = "ord") %>%
  ## Fill the NAs for the unknown families in these orders
  dplyr::mutate(probem = ifelse(is.na(probem) & is.na(family),
                                probem_ord, probem)) %>%
  ## Rename column
  dplyr::select(-probem_ord)

  
# Pull out known families that are still NA
bwg_probem %>% 
  dplyr::filter(is.na(probem)) %>% 
  dplyr::select(family) %>% 
  dplyr::distinct()

# Get references for these
## Odonata, hard to predict, not tropical values, use temperate upper value
## because eggs not laid in large batches. Survival likely higher
### Lawton, J.H. A population study on larvae of the damselfly Pyrrhosoma 
### nymphula (Sulzer) (Odonata: Zygoptera). Hydrobiologia 36, 33–52 (1970). 
### https://doi.org/10.1007/BF00751281
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(!is.na(ord) & ord == "Odonata",
                                0.3, probem))


## These families are semi-aquatic/terrestrial, hard to say so use 0
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(!is.na(family) & family %in% 
                                  c("Carabidae", "Curculionidae", 
                                    "Elateridae", "Enicocephalidae",
                                    "Tineidae", "Scutelleridae"),
                                0, probem))



## Calamoceratidae
### No data on family
### Fan, Jijun, Tan, Lin, Yang, Yishuai, Liu, Wei, Tan, Jihui, Li, Yong, Miao, 
### Youcheng, Hu, Qiulong and Zhou, Linghong. "Temperature-dependent development, 
### predation, and life table of Sphaerophoria macrogaster (Thomson) 
### (Diptera: Syrphidae) feeding on Myzus persicae (Sulzer) (Homoptera: Aphididae)" 
### Open Life Sciences, vol. 20, no. 1, 2025, pp. 20251183. 
### https://doi.org/10.1515/biol-2025-1183
### Here the only study is laboratory,gives about 60% percent survival. 
### Nothing else to estimate so say that field survival is half that (predators, 
### environmental fluctuations...)
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(!is.na(ord) & ord == "Trichoptera",
                                0.29, probem))

## Use Coleoptera because no reference
### Lampyridae artificial rearing studies
### https://agris.fao.org/search/en/records/6472525d53aa8c8963062412?utm_source=chatgpt.com
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(!is.na(family) & family %in% 
                                  c("Lampyridae", "Ptilodactylidae"),
                                ### only one value for coleopterans
                                ol_probem %>% 
                                  dplyr::filter(species == "Cyphon") %>% 
                                  dplyr::select(probem) %>% 
                                  dplyr::pull(), probem))

## Use Diptera because no reference
### NB Ephydridae has one study but in lab and survical very high (90%+)
### Ugine TA, Sanderson JP, Wraight SP. Developmental times and life tables for 
### shore flies, Scatella tenuicosta (Diptera: Ephydridae), at three temperatures. 
### Environ Entomol. 2007 Oct;36(5):989-97. doi: 10.1603/0046-225x(2007)36[989:dtaltf]2.0.
### Stratomyidae only controlled studies on Hermetia illucens
### Syrphidae only studies in controlled lab environments
### amdem, M. M., & Otomo, P. V. (2023). Developmental performance of Eristalis 
### tenax larvae (Diptera: Syrphidae): Influence of growth media and yeast addition 
### during captive rearing. Journal of Experimental Zoology Part A: Ecological and 
### Integrative Physiology, 339, 503–513. https://doi.org/10.1002/jez.2696
### aff. Drosophilidae because ID uncertain
bwg_probem <- 
  bwg_probem %>% 
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
                                ### Extract average of all Dipterans
                                ol_probem %>% 
                                  dplyr::filter(species != "Cyphon") %>% 
                                  dplyr::select(probem) %>% 
                                  dplyr::summarise(dplyr::across(probem,
                                                                 ~mean(.x))) %>% 
                                  dplyr::pull(), probem))

# Now what is left is Hemiptera, Ephemeroptera, Lepidoptera, all unided, put 0
bwg_probem <- 
  bwg_probem %>% 
  dplyr::mutate(probem = ifelse(is.na(probem),
                                0, probem))

# save for email
readr::write_csv(bwg_probem,
                 here::here("bwg_probemergence.csv"))




  