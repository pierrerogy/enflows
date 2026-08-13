# Functions
# Author: Pierre Rogy

# Libraries
library(tidyverse)
library(curl)
library(here)
library(cheddar)
library(viridis)
library(brms)
library(latex2exp)
library(ggrepel)

# Negate operator
`%notin%` <- 
  Negate(`%in%`)

# Make matrices binary ----------------------------------------------------
make_binary_matrices <- function(matrix_list,
                                 cutoff){
# Make empty list to store the new ones in
ret <- 
  list()
# Loop with cutoff
for (i in 1:length(matrix_list)){
  mat <- 
    as.matrix(matrix_list[[i]])
  ## All positive values above cutoff are 1
  mat[mat >= cutoff] <- 
    1
  ## Rest is 0
  mat[mat < cutoff] <- 
    0 
  ## Save to list of matrices
  ret[[i]] <- 
    mat
}
# Keep names
names(ret) <- 
  names(matrix_list)

# Return
return(ret)

}


# Make food webs -----------------------------------------------------------
make_food_webs <- function(binary_matrices, biomass_data){
  ## Mostly written by Eoin O'Gorman
  # This line is really important or character text will be seen as factors, 
  # creating problems in Cheddar
  options(stringsAsFactors = FALSE)   
  
  # Create tibble that will store everything
  ret <- 
    tibble::tibble()
  
  # Loop to compute all foodwebs
  for(i in 1:length(binary_matrices)){
    ## Select relevant matrix
    site_mat <- 
      binary_matrices[[i]]
    ## Get site name
    sitename <- 
      names(site_matrices)[i]
    
    # Make vector of bromeliads in site
    brom_list <- 
      biomass_data %>%
      dplyr::filter(site_code == sitename) %>% 
      dplyr::select(bromeliad_id) %>% 
      unique() %>% 
      dplyr::pull()
    
    ## Make food web for each bromeliad at that site 
    for(j in brom_list){
      ## Filter biomass data by bromeliad
      brom_biomass <- 
        biomass_data %>% 
        dplyr::filter(bromeliad_id == j)
      ## Get species list
      speclist <- 
        brom_biomass %>% 
        dplyr::select(bwg_name) %>% 
        unique() %>% 
        dplyr::pull()
      
      ## Subset site matrix by what is inside the bromeliad
      brom_mat <- 
        site_mat[rownames(site_mat) %in% speclist,
                 colnames(site_mat) %in% speclist]
      
      ## Convert to trophic link
      trophic_links <- 
        data.frame(Matrix.to.list(brom_mat)) %>% 
        dplyr::rename(consumer = X1,
                      resource = X2) %>% 
        dplyr::filter(consumer != "0") %>% 
        dplyr::mutate(site_code = sitename,
                      bromeliad_id = j)
      
      ## Now if we only have the two zooplankton consumers make empty
      if(length(unique(trophic_links$consumer)) == 2){
        ### make empty
        trophic_links <- 
          tibble::tibble()
      }
      
      ## Add to tibble
      ret <- 
        ret %>% 
        bind_rows(trophic_links)
    }
    
    
  }
  
  # Return
  return(ret)
  
}

# Take a food web in matrix format and converts it to list format ---------
## By Andrew MacDonald, modified by me
## (plotting of food webs)
Matrix.to.list <- function(web.matrix, predator.first = TRUE){
  if(length(dimnames(web.matrix)[[1]]) == length(web.matrix[1,]))
    species.names <- dimnames(web.matrix)[[1]] else
      species.names <- 1:length(web.matrix[,1])
    web.list <- matrix(0, sum(web.matrix), 2)
    counter <- 1
    for(i in 1:length(web.matrix[,1]))
      for(j in 1:length(web.matrix[,1]))
        if(web.matrix[j,i]==1){   ##Note: swapped i and j to fix problem with order of predator/prey columns
          web.list[counter,] <- c(species.names[i],species.names[j])
          counter <- counter + 1
        }
    if(!predator.first)
      web.list <- cbind(web.list[,2], web.list[,1])
    web.list
}



# Load matrices from a folder ---------------------------------------------
load_matrices <- function(folder, ending){
  
  # Load filenames
  filenames <- 
    dir(folder) 
  
  # Make empty list to store matrices
  ret <- 
    vector("list",
           length = length(filenames))
  
  # Name each element
  names(ret) <- 
    filenames %>% 
    stringr::str_remove(ending)
  
  # Extract matrices
  for (k in 1: length(filenames)){
    ret[[k]] <- 
      read.table(paste0(folder, "/",
                        filenames[[k]]), 
                 sep = ',',
                 header = T,
                 row.names = 1
      ) 
  }
  
  # Return
  return(ret)
  
}


# Convert matrices to long format ------------------------------------------------
make_mat_long <- function(matrix, unit){
  
  # Convert
  ret <-
    matrix %>% 
    ## Convert to tibble
    tibble::tibble() %>% 
    ## Create resource column
    dplyr::mutate(resource = colnames(.)) %>% 
    ## Make data long format
    tidyr::pivot_longer(cols = where(is.numeric),
                        values_to = unit,
                        names_to = "consumer")
  
  # Return
  return(ret)

}
# Get CP3 values ----------------------------------------------------------
match_CP3 <- function(ord, subord, family, genus, cp3_search){
  # Return object
  ret <- numeric(0)
  
  # Cycle through taxa in increasing order
    for(target in c(genus, family, subord, ord)){
      
      ### Filter and extract taxon in cp3_search
      ret <- 
        cp3_search %>% 
        dplyr::filter(taxon == target) %>% 
        dplyr::select(cp3_val) %>% 
        dplyr::pull()
    
      ## Break if we got the value we wanted
      if(length(ret) > 0 ){
        break
        
      }
      
    }
    
  # Return
  return(ret)
  
  }
  


# Plot PCA outputs --------------------------------------------------------
plz_plot_my_pca <- function(PCA, axes, cor_min = 0.5, type, dats){
  

  
  # Get proportion explained for each axis and convert to %
  per_inertia <- 
    (PCA$eigenvalues[axes,2]) * 100
  
  if(type == "points"){
    # Make a taxon vector to make color points
    taxa_name <-
      ## Get fw data
      dats %>% 
      ## Make new column
      dplyr::mutate(taxa_name = ord) %>% 
      ## Fill in the blanks
      dplyr::mutate(taxa_name = ifelse(is.na(taxa_name) & !is.na(subclass),
                                       subclass, 
                                       ifelse(is.na(taxa_name) & !is.na(class),
                                              class, 
                                              ifelse(is.na(taxa_name) & !is.na(phylum),
                                                     phylum, 
                                                     taxa_name)))) %>% 
      ## Get main Diptera families in
      dplyr::mutate(taxa_name = ifelse(family == "Chironomidae"& !is.na(family),
                                       "Chironomidae",
                                       ifelse(family == "Culicidae" & !is.na(family),
                                              "Culicidae",
                                              taxa_name))) %>%
      ## Get other Diptera families in, and other taxonomic groups
      dplyr::mutate(taxa_name = ifelse(taxa_name == "Diptera",
                                       "Other_Diptera", 
                                       ifelse(taxa_name == "Trombidiformes", 
                                              "Acari",
                                              ifelse(taxa_name == "Rhynchobdellida", 
                                                     "Hirudinea",
                                                     ifelse(taxa_name == "Turbellaria", 
                                                            "Platyhelminthes",
                                                            ifelse(taxa_name %in% c("Harpacticoida", "Cyclopoida"),
                                                                   "Copepoda",
                                                                   ifelse(taxa_name == "Haplotaxida", 
                                                                          "Oligochaeta",
                                                                          ifelse(taxa_name == "Nematomorpha", 
                                                                                 "Nematoda", taxa_name)))))))) 
    
    
    # Select axes 
    ind_scores <-
      PCA$individuals[,axes] %>%
      tibble::as_tibble() %>% 
      dplyr::rename("first" = 1,
                    "second" = 2) %>% 
      ## add taxa name
      dplyr::bind_cols(taxa_name %>% 
                         dplyr::select(taxa_name)) %>% 
      ## Add colours
      dplyr::mutate(col = ifelse(taxa_name == "Chironomidae",
                                 "black",
                                 ifelse(taxa_name == "Culicidae",
                                        "grey70",
                                        ifelse(taxa_name == "Other_Diptera",
                                               "grey30",
                                               ifelse(taxa_name %in% c("Opisthophora", "Oligochaeta",
                                                                       "Hirudinea","Platyhelminthes"),
                                                      "skyblue",
                                                      ifelse(taxa_name == "Coleoptera",
                                                             "firebrick",
                                                             ifelse(taxa_name == "Odonata",
                                                                    "mediumorchid4",
                                                               ifelse(taxa_name %in% c("Lepidoptera", "Trichoptera"),
                                                                      "darkgreen",
                                                                      ifelse(taxa_name %in% c("Hemiptera", "Acari", "Panpulmonata"),
                                                                             "yellowgreen",
                                                                             ifelse(taxa_name %in% c("Diplostraca", "Copepoda", "Podocopida"),
                                                                                    "steelblue4", NA))))))))),
                    pch = ifelse(taxa_name %in% c("Opisthophora", "Oligochaeta",
                                                  "Hirudinea","Platyhelminthes"),
                                 17,
                                 ifelse(taxa_name %in% c("Chironomidae", "Culicidae", "Other_Diptera"),
                                        15,
                                        ifelse(taxa_name %in% c("Odonata",  "Coleoptera", "Lepidoptera", "Trichoptera"),
                                               19,
                                               ifelse(taxa_name %in% c("Hemiptera", "Acari", "Panpulmonata"),
                                                      4, 
                                                      ifelse(taxa_name %in% c("Diplostraca", "Copepoda", "Podocopida"),
                                                             18,
                                                             NA))))))
    
    # Get limits for plot
    axis_limits <- 
      c(min(ind_scores[,1:2]), 
        max(ind_scores[,1:2]))
    
    # Make plot
    ret <-
    ## Set plot
    ggplot(aes(x = first,
               y = second, 
               colour = col,
               shape = pch),
           data = ind_scores)+
      ## Styling
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.title=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA)) +
      ## Add points
      geom_point(size = 3) +
      ## Colour is exactly what we are giving it
      scale_shape_identity(guide = "legend",
                           labels = c("Diptera", 
                                      "Oligochaetea, Nematoda,\nHirudinea", 
                                      "Other insects", 
                                      "Terrestrials", 
                                      "Crustaceans"),
                           breaks = c(15, 17, 19, 4, 18)) +
      scale_colour_identity(guide = "legend",
                            labels = c("Chironomidae",
                                       "Culicidae",
                                       "Other Diptera",
                                       "Oligochaetea, Nematoda,\nHirudinea",
                                       "Coleoptera",
                                       "Odonata",
                                       "Other insects",
                                       "Terrestrials",
                                       "Crustaceans",
                                       ""),
                            breaks = c("black", 
                                       "grey70", 
                                       "grey30", 
                                       "skyblue",
                                       "firebrick", 
                                       "mediumorchid4",
                                       "darkgreen",  
                                       "yellowgreen",
                                       "steelblue4",  
                                       NA))}
  
  if(type == "arrows"){
    # Get scores of variables
    var_scores <- 
      PCA$variables[,axes]
    
    # Keep variables above correlation cutoff
    var_scores <-
      var_scores[apply(abs(var_scores) >= cor_min, 
                       1, sum ) >0 ,]  
    
    # Define circle
    circle <- 
      tibble::tibble(x = cos(seq(0, 2*pi, length = 200)),
                     y = sin(seq(0, 2*pi, length = 200)))
  
    # Plot!
    ret <-
      ## Set circle
      ggplot(aes(x = x,
                 y = y),
             data = circle) +
      ## Add circle line
      geom_path() +
      ## Styling
      theme(legend.justification = NULL,
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line =  element_blank()) +
      ## Add all arrows to circle plot
      geom_segment(aes(x = 0, 
                       y = 0, 
                       xend = var_scores[,1]*1, 
                       yend = var_scores[,2]*1),
                   data = var_scores,
                   arrow = arrow()) +
      ## add arrow labels based on row names
      ggrepel::geom_text_repel(aes(x = (var_scores[,1]*1.1),
                                   y = (var_scores[,2]*1.1),
                                   label = rownames(var_scores)),
                               data = var_scores,
                               size = 3)}
  
  # Last customisation of plot
  ret <- 
    ret +
    ## Dotted x axis
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "gray") +
    ## Dotted x axis
    geom_vline(xintercept = 0,
               linetype = "dashed",
               color = "gray") +
    ## X axis label based on % explained
    xlab(paste0("PC ", 
                axes[1], 
                " (", round(per_inertia[1], 2), 
                "%)")) +
    ## Y axis label based on % explained
    ylab(paste0("PC ", 
                axes[2], 
                " (", 
                round(per_inertia[2], 2), 
                "%)")) 
  # Return
  return(ret)
}

# Food web functions ------------------------------------------------------
## Trophic level by Benoît Gauzens
TL <- function(LL)
{
  LL <- t(LL)
  nn <- rowSums(LL); nn[nn==0] <- 1
  ww <- diag(1/nn)
  L1 <- ww %*% LL
  L2 <- L1 - diag(rep(1,length(nn)))
  b <- -1*rep(1,length(nn))
  Tro.lev <- solve(L2) %*% b
  return(Tro.lev)
}



# Get vector of x-axis positions ------------------------------------------
get_TL_pos <- function(dats){
  
  # Tibble on which we will add the positions
  ret <- 
    tibble::tibble()
  
   # We are bounded between -1 and 1
  # For each trophic level
  for(i in unique(dats$TL)){
    ## Make new column with equal space between all species
    temp <- 
      dats %>% 
      dplyr::filter(TL == i) %>% 
      dplyr::mutate(pos = seq(from = -1, 
                              to = 1, 
                              by = 2 / (nrow(.) - 1)))
    ## Don't know how to mathematise this..
    if(nrow(temp) == 2){
      temp <- 
        temp %>% 
        dplyr::mutate(pos = c(-0.5, 0.5))
    }
    
    ## Add it to return
    ret <- 
      ret %>% 
      dplyr::bind_rows(temp)

  }
  
  # Clean and return
  dats %>% 
    dplyr::left_join(ret,
                     by = c("species", "TL")) %>% 
    dplyr::select(pos) %>% 
    dplyr::pull() %>% 
    return()

}

# Plot food web -----------------------------------------------------------
plot_food_web <- function(bromtoplot,
                          sitetoplot,
                          fluxes_long_format,
                          clean_biomass_data,
                          binary_matrices){
  
  # Define basal species
  basal <- 
    c("algae", "FPOM", "CPOM", "microzooplankton",
    "mesozooplankton", "bromeliad.tissue", "microorganisms")
  
  # # Left aside for now
  # # Define functional group shape coding
  # shapes <- 
  #   c("resource" = NA, "scraper"= 15 , "gatherer" = 16, "filter.feeder" = 17, 
  #     "shredder" = 18, "predator" = 4 , "unknown" = 13)
  
  # Filter bromeliad in web
  web <- 
    fluxes_long_format %>% 
    dplyr::filter(site_code == sitetoplot & bromeliad_id == bromtoplot) %>% 
    ## Remove cannibalism
    dplyr::filter(consumer != resource)
    
  
  # Get list of species in the web
  speclist <- 
    unique(c(web$consumer,
             web$resource))
  
  # Get site interaction matrix
  site_mat <- 
    binary_matrices[[which(names(binary_matrices) == sitetoplot)]]
  
  # Get trophic levels within this bromeliad
  levels <- 
    ## Subset site matrix by what is inside the bromeliad
    site_mat[rownames(site_mat) %in% speclist,
             colnames(site_mat) %in% speclist] %>% 
    ## Use function to get trophic levels
    TL() %>% 
    tibble::as_tibble(rownames = "species",
                      .name_repair = "minimal") %>% 
    dplyr::rename(TL = 2) %>% 
    ## Give random x axis position for organisms
    ## MODIFY WHEN I GET TIME
    dplyr::mutate(pos = runif(n = nrow(.), min = -1, max = 1)) %>% 
    ## Give specific ones for resources in resource column
    dplyr::mutate(pos = ifelse(species == "CPOM",
                               -1, ifelse(species == "FPOM",
                                          -0.6, ifelse(species == "microorganisms",
                                                       -0.3, ifelse(species == "microzooplankton",
                                                                    0, ifelse(species == "mesozooplankton",
                                                                              0.3, ifelse(species == "algae",
                                                                                          0.6, ifelse(species == "bromeliad.tissue",
                                                                                                      1, pos)))))))) %>% 
    ## Add biomass
    dplyr::left_join(clean_biomass_data %>% 
                       dplyr::filter(site_code == sitetoplot & bromeliad_id == bromtoplot) %>% 
                       dplyr::select(biomass = tot_biomass_g,
                                     species = bwg_name),
                     by = "species") %>% 
    ## Add label, done inefficiently
    dplyr::mutate(
      ## Remove first four characters of string
      label = stringr::str_remove(string = species,
                                  pattern = "^[a-zA-Z]{4}"),
      ## Remove dot and 1-3 numerical characters after it
      label = stringr::str_remove(string = label,
                                  pattern = "\\.[0-9]{1,3}$"),
      ## To prevent warnings from happening, replace empty strings
      label = ifelse(label == "",
                     "blablabla", label)) %>% 
    ## Now remove this string from species name 
    dplyr::mutate(label = stringr::str_remove(string = species,
                                              pattern = label)) %>% 
    ## Custom labels for some resources
    dplyr::mutate(label = ifelse(species == "algae",
                                 "algae", 
                                 ifelse(species == "FPOM", 
                                        "FPOM", 
                                        ifelse(species == "CPOM",
                                               "CPOM",
                                               ifelse(species == "microorganisms",
                                                      "bact", label))))) 
  
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
  
  # # Last step finally, adding functional group to web
  # # GIVE UP FOR NOW
  # levels <- 
  #   levels %>% 
  #   ## Add functional group
  #   dplyr::left_join(web %>% 
  #                      dplyr::select(consumer, 
  #                                    fg = cons_functional_group) %>% 
  #                      dplyr::distinct() %>% 
  #                      dplyr::rename(species = consumer),
  #                    by = "species") %>% 
  #   ## add resource functional group, and add exception for unknown ones
  #   dplyr::mutate(fg = ifelse(species %in% basal & is.na(fg),
  #                             "resource", 
  #                             ifelse(species %notin% basal & is.na(fg),
  #                                    "unknown", fg))) %>% 
  #   ## finally make piercers and engulfers as a single predator group
  #   dplyr::mutate(fg = ifelse(fg %in% c("engulfer", "piercer"),
  #                             "predator", fg)) %>% 
  #   ## aaaand add shape parameter
  #   dplyr::left_join(shapes,
  #                    by = "fg")
  
  
  # Plot
  ## Base plot on levels tibble
  ret <- 
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
      # ## Add point shape based on functional group kind of gave up because 
      # ## adds to much, kinda ugly
      # geom_point(aes(shape = fg), 
      #            size = 7) +
      # scale_shape_manual(values = shapes,
      #                    name = "Functional \ngroup")
      #geom_text(data = colvec_x, aes(x = x, y = y, label = code), colour = 'black', size = 2) +
      ## Species labels
      geom_label(aes(label = label),
                 color = "white",
                 fill = "black",
                 size = 3) +
      ## Styling
      theme(legend.justification = NULL,
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black"))
  
  # Return plot
  return(ret)
}




# Labels for plots --------------------------------------------------------
full_label <- function(var){
  
  
  # First set of hypotheses
  if(var == "flux_raw")
    ret <- 
      expression(paste("Energy flux entering a \nconsumer node (log J.y"^"-1"*")"))
  if(var == "Axis.1")
    ret <- 
      "Length of life cycle (PC1 score)"
  if(var == "Axis.2")
    ret <- 
      "Pelagicity (PC2 score)"
  if(var == "Axis.3")
    ret <- 
      "Morphological defense (PC3 score)"
  if(var == "Axis.4")
    ret <- 
      "Simplicity of life cycle (PC4 score)"
  
  # Second set of hypotheses
  if(var == "flux_raw_pair")
    ret <- 
      expression(paste("Energy flux between a \nconsumer-resource pair (log J.y"^"-1"*")"))
  if(var == "dist_axis1")
    ret <- 
      "Difference in length of life cycle \n(PC1 scores)"
  if(var == "dist_axis2")
    ret <- 
      "Difference in pelagicity \n(PC2 scores)"
  if(var == "dist_axis3")
    ret <- 
      "Difference in morphological defense \n(PC3 scores)"
  if(var == "dist_axis4")
    ret <- 
      "Difference in simplicity of life cycle \n(PC4 scores)"
  if(var == "func_dist")
    ret <- 
      "Euclidean distance in fuctional \nhypervolume"
  
  
  # Third set of hypotheses
  if(var == "gini")
    ret <- 
      "Gini index of energy fluxes"
  if(var == "totalflow_Jy")
    ret <- 
      expression(paste("Total energy flow (log J.y"^"-1"*")"))
  if(var == "predflow_Jy")
    ret <- 
      expression(paste("Energy flow to \npredators (log J.y"^"-1"*")"))
  if(var == "prop_hull")
    ret <- 
      expression(paste("Proportion of overall functional \nhull covered by food web"))
  if(var == "avg_func_dist")
    ret <- 
      "Average Euclidean distance \nin fuctional hypervolume"

  # Return
  return(ret)
}

# Basic line plot ---------------------------------------------------------
basic_line <- function(model, dats, x_var, y_var, index = 1, scale = NA){
  # Get effect
  model_effect <- 
    brms::conditional_effects(model,
                              method = "fitted")[[index]]

  # x axis not logged
  if(scale == "x_scale"){
    ret <- 
      ggplot(data = model_effect,
           aes(x = .data[[x_var]],
               y = estimate__)) +
      ## Add raw data points
      geom_jitter(data = dats,
                 aes(x = .data[[x_var]],
                     y = log(.data[[y_var]])),
                 size = 0.5,
                 width = 0.001,
                 alpha = 0.4)} 
  
  # y axis not logged
  if(scale == "y_original"){
    ret <- 
      ggplot(data = model_effect,
             aes(x = log(.data[[x_var]]),
                 y = estimate__)) +
      ## Add raw data points
      geom_jitter(data = dats,
                  aes(x = log(.data[[x_var]]),
                      y = .data[[y_var]]),
                  size = 0.5,
                  width = 0.001,
                  alpha = 0.4)} 
  
  # Everything was logged
  if(scale == "none"){ 
    ret <- 
      ggplot(data = model_effect,
             aes(x = log(.data[[x_var]]),
                 y = estimate__)) +
      ## Add raw data points
      geom_jitter(data = dats,
                  aes(x = log(.data[[x_var]]),
                      y = log(.data[[y_var]])),
                  size = 0.5,
                  width = 0.001,
                  alpha = 0.4)}
  
  # Nothing was logged
  if(scale == "both_original"){ 
    ret <- 
      ggplot(data = model_effect,
             aes(x = .data[[x_var]],
                 y = estimate__)) +
      ## Add raw data points
      geom_jitter(data = dats,
                  aes(x = .data[[x_var]],
                      y = .data[[y_var]]),
                  size = 0.5,
                  width = 0.001,
                  alpha = 0.4)}
  
  # Common elements to plots
  ret <- 
    ret +
    ## Lines
    geom_line(linewidth = 1.2,
              colour = "steelblue") +
    ## Confidence interval
    geom_ribbon(aes(ymin = lower__,
                    ymax = upper__),
                alpha = 0.5,
                colour = NA,
                fill = "steelblue") +
    ## Styling
    xlab(full_label(x_var)) +
    ylab(full_label(y_var)) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white"))                
  
  # Return plot
  return(ret)
  
  
  
  
}

# Get title for significant relationships ---------------------------------
title_sig <- function(model, index, interaction = FALSE){
  # Helper - pull the pd and % in ROPE for one parameter, format the label,
  # and flag whether that parameter clears the significance threshold
  # (pd greater than 99.9% and less than 5% of the posterior inside the ROPE)
  summarise_term <- function(param_index){
    ## Probability of direction, as a percentage
    pd <-
      round(bayestestR::pd(model)$pd[param_index],
            digits = 3) * 100
    ## Percentage of the posterior falling within the ROPE
    rope <-
      round(bayestestR::rope(model)$ROPE_Percentage[param_index],
            digits = 3) * 100
    ## Return the formatted text and the significance flag together
    list(text = paste0("pd = ",
                       round(pd,
                             digits = 3),
                       "%; ",
                       "% rope = ",
                       round(rope,
                             digits = 3)),
         sig = pd > 99.9 & rope < 5)
  }

  # In cases where there is no interaction to plot --------------------------
  if(isFALSE(interaction)){
    ## Summarise the single focal term (e.g. the quadratic term at index 3)
    term <-
      summarise_term(index)
    ## Keep the label, and set bold face only if the term is significant
    ret <-
      list(label = term$text,
           fontface = ifelse(term$sig,
                             "bold",
                             "plain"))

  # In cases where we have to plot the interaction --------------------------
  } else {
    ## Initiate the label string and a running significance flag
    label <-
      ""
    sig <-
      FALSE

    ## Two reps of extracting coefficients (resource higher, then lower)
    for(i in c(2, 4)){
      ## Add beginning of string depending on rep
      if(i == 2)
        label <-
          paste0(label, "R. higher: ")
      if(i == 4)
        label <-
          paste0(label, ";\nR. lower: ")

      ## Summarise this arm of the interaction
      term <-
        summarise_term(i)
      ## Append its text, and flag significance if either arm is significant
      label <-
        paste0(label, term$text)
      sig <-
        sig | term$sig}

    ## Keep the label, and set bold face if either arm is significant
    ret <-
      list(label = label,
           fontface = ifelse(sig,
                            "bold",
                            "plain"))}

  ## Return
  return(ret)
}

# Fit the whole suite of models for a single site -------------------------
fit_site_suite <- function(site, percons, perpair, perbrom){
  # Print site name as a sanity check while the suite runs
  print(site)

  # Subset the per-consumer data to the focal site
  dats_percons <-
    percons %>%
    dplyr::filter(site_code == site)
  # Subset the per-pair data to the focal site
  dats_perpair <-
    perpair %>%
    dplyr::filter(site_code == site)
  # Subset the per-bromeliad data to the focal site
  dats_perbrom <-
    perbrom %>%
    dplyr::filter(site_code == site)

  # Return the whole suite as a named list
  list(
    # Hypothesis 1 - raw flow entering a consumer node ~ each PC axis -----
    ## Axis 1
    consflow_model_axis1 <- 
      brms::brm(log(flux_raw + 0.001) ~
                  scale(Axis.1,
                        scale = FALSE, center = TRUE) +
                  I(scale(Axis.1, scale = FALSE, center = TRUE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.99,
                               max_treedepth = 15),
                data = dats_percons),
    ## Axis 2
    consflow_model_axis2 <- 
      brms::brm(log(flux_raw + 0.001) ~
                  scale(Axis.2,
                        scale = FALSE, center = TRUE) +
                  I(scale(Axis.2, scale = FALSE, center = TRUE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.99,
                               max_treedepth = 15),
                data = dats_percons),
    ## Axis 3
    consflow_model_axis3 <- 
      brms::brm(log(flux_raw + 0.001) ~
                  scale(Axis.3,
                        scale = FALSE, center = TRUE) +
                  I(scale(Axis.3, scale = FALSE, center = TRUE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.99,
                               max_treedepth = 15),
                data = dats_percons),
    ## Axis 4
    consflow_model_axis4 <- 
      brms::brm(log(flux_raw + 0.001) ~
                  scale(Axis.4,
                        scale = FALSE, center = TRUE) +
                  I(scale(Axis.4, scale = FALSE, center = TRUE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.99,
                               max_treedepth = 15),
                data = dats_percons),

    # Hypothesis 2 - raw flow in a pair ~ functional distance ------------
    ## Overall functional distance
    pairflow_model_overall <- 
      brms::brm(log(flux_raw_pair + 0.001) ~
                  scale(func_dist,
                        center = TRUE, scale = FALSE) +
                  I(scale(func_dist, center = TRUE, scale = FALSE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perpair),
    ## Distance along Axis 1
    pairflow_model_axis1 <- 
      brms::brm(log(flux_raw_pair + 0.001) ~
                  scale(dist_axis1,
                        center = TRUE, scale = FALSE) +
                  I(scale(dist_axis1, center = TRUE, scale = FALSE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perpair),
    ## Distance along Axis 2
    pairflow_model_axis2 <- 
      brms::brm(log(flux_raw_pair + 0.001) ~
                  scale(dist_axis2,
                        center = TRUE, scale = FALSE) +
                  I(scale(dist_axis2, center = TRUE, scale = FALSE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perpair),
    ## Distance along Axis 3
    pairflow_model_axis3 <- 
      brms::brm(log(flux_raw_pair + 0.001) ~
                  scale(dist_axis3,
                        center = TRUE, scale = FALSE) +
                  I(scale(dist_axis3, center = TRUE, scale = FALSE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perpair),
    ## Distance along Axis 4
    pairflow_model_axis4 <- 
      brms::brm(log(flux_raw_pair + 0.001) ~
                  scale(dist_axis4,
                        center = TRUE, scale = FALSE) +
                  I(scale(dist_axis4, center = TRUE, scale = FALSE)^2) +
                  (1|bromeliad_id),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perpair),

    # Hypothesis 3 - whole-web flow ~ functional volume and distance -----
    ## Total flow against hull volume
    totalflow_model_vol <- 
      brms::brm(log(totalflow_Jy) ~
                  scale(log(prop_hull),
                        center = TRUE, scale = FALSE) +
                  I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perbrom),
    ## Total flow against average functional distance
    totalflow_model_dist <-
      brms::brm(log(totalflow_Jy) ~
                  scale(avg_func_dist,
                        center = TRUE, scale = FALSE) +
                  I(scale(avg_func_dist, center = TRUE, scale = FALSE)^2),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perbrom),
    ## Flow to predators against hull volume
    predflow_model_vol <-
      brms::brm(log(predflow_Jy + 0.01) ~
                  scale(log(prop_hull),
                        center = TRUE, scale = FALSE) +
                  I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perbrom),
    ## Flow to predators against average functional distance
    predflow_model_dist <- 
      brms::brm(log(predflow_Jy + 0.01) ~
                  scale(avg_func_dist,
                        center = TRUE, scale = FALSE) +
                  I(scale(avg_func_dist, center = TRUE, scale = FALSE)^2),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perbrom),
    ## Gini index of fluxes against hull volume
    gini_model_vol <- 
      brms::brm(gini ~
                  scale(log(prop_hull),
                        center = TRUE, scale = FALSE) +
                  I(scale(log(prop_hull), center = TRUE, scale = FALSE)^2),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perbrom),
    ## Gini index of fluxes against average functional distance
    gini_model_dist <- 
      brms::brm(gini ~
                  scale(log(avg_func_dist),
                        center = TRUE, scale = FALSE) +
                  I(scale(log(avg_func_dist), center = TRUE, scale = FALSE)^2),
                iter = 20000,
                family = gaussian(),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 15),
                data = dats_perbrom))
}


