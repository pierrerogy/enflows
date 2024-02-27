# Functions
# Author: Pierre Rogy

# Libraries
library(tidyverse)
library(curl)
library(here)
library(raster)
library(sp)
library(cheddar)
library(viridis)
library(brms)

# Download data from chelsa -----------------------------------------------
get_chelsa <- function(txt_file, output_folder){
  # First read list of files from txt file
  files <- 
    readr::read_table(txt_file,
                      col_names = F) %>% 
    dplyr::select(X1) %>% 
    dplyr::pull()
  
  # Download files
  for(i in files){
    # Increase download limit
    options(timeout = max(NA, getOption("timeout")))
    curl::curl_download(url = i,
                        mode = "wb",
                        destfile = paste0(output_folder,
                                    "/",
                                    ## Here just keep text after the last / as file name
                                    stringr::str_match(string = i,  
                                                       pattern = "([^/]+$)")[1]))
  }
  }

# Load the raster files and extract the information of the raster -------------------------------------------------------------------------
get_bio_data <- function(raster_files, sites_coords){
  
  # Rename supplied data
  ret <- 
    sites_coords
  
  # Convert coords to spatial points
  coords <- 
    sp::SpatialPoints(ret[,2:3])
  
  # Make vector of variable names
  var_names <- 
    c("site", "lon", "lat")
  
  # For each raster file
  for(i in raster_files){
  
  # Load raster file
  ras <- 
    raster::raster(i)
  
  # Get name of variable
  ras_name <- 
    names(ras) %>% 
    stringr::str_split(pattern = "_") %>% 
    purrr::flatten_chr() %>% 
    dplyr::nth(2)
  
  # Add name to list of names
  var_names <- 
    c(var_names,
      ras_name)
  
  # Extract values
  ret <-
    ret %>% 
    dplyr::bind_cols(raster::extract(ras,
                                     coords))
  
  }
  
  # Give correct column names
  colnames(ret) <- 
    var_names
  
  
  # Return
  return(ret)
}


# Add environmental variables to dataframe --------------------------------
add_env <- function(dats, raster_folder){
  # Get list of raster files
  raster_files <- 
    paste0(raster_folder, 
           "/",
           list.files(raster_folder))
  
  # Extract site info
  sites_coords <- 
    dats %>% 
    dplyr::select(site, lon, lat) %>% 
    dplyr::distinct()
  
  # Extract chelsa values and attach to coordinates
  env_values <- 
    get_bio_data(raster_files = raster_files, 
                 sites_coords = sites_coords)
  
  # Stick to original data
  ret <- 
    dats %>% 
    dplyr::left_join(env_values,
                     by = c("site", "lat", "lon"))
  
  # Return
  return(ret)
  
}



# not in ------------------------------------------------------------------
'%notin%' <- 
  Negate('%in%')



# Get equivalent models from loo output-------------------------------------------------------------------------
get_eq_models <- function(loo_output){
  ret <- 
    ## Extract diff table values
    tibble::tibble(model = attr(loo_output, "dimnames")[[1]],
                   elpd_diff = loo_output[,1],
                   se_diff = loo_output[,2]) %>% 
    ## Get ratio of elpd/se to get equivalent models
    dplyr::mutate(ratio = abs(elpd_diff/se_diff)) %>% 
    ## Keep best and equivalent models
    dplyr::filter(ratio < 2 | is.nan(ratio))
  
  # Return
  return(ret)
}

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
make_food_webs <- function(taxon_names, basal_species, binary_matrices, biomass_data){
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
      biomass_clean %>%
      dplyr::filter(site_code == sitename) %>% 
      dplyr::select(bromeliad_id) %>% 
      unique() %>% 
      dplyr::pull()
    
    ## Make food web for each bromeliad at that site 
    for(j in brom_list){
      ## Filter biomass data by bromeliad
      brom_biomass <- 
        biomass_clean %>% 
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
        ### Print to check
        print(unique(trophic_links$consumer))
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


# Convert matrices to long ------------------------------------------------
make_mat_long <- function(matrix, unit){
  
  # Convert
  ret <-
    matrix %>% 
    ## Convert to tibble
    tibble() %>% 
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
    dplyr::filter(site == sitetoplot & bromeliad_id == bromtoplot) %>% 
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
    tibble::as_tibble(rownames = "species") %>% 
    dplyr::rename(TL = V1) %>% 
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
full_label <- function(var, scale = "none"){
  
  if(var == "cons_P")
    ret <- 
      "P content of consumer"
  
  if(var == "flux_prop")
    ret <- 
      "Prop. of total energy flux \nentering consumer node"
  
  if(var == "flux_raw")
    ret <- 
      expression(paste("Energy flux entering a \nconsumer node (log J.y"^"-1"*")"))
  
  if(var == "flux_prop_pair")
    ret <- 
      "Prop. of total energy flux \nbetween a consumer-resource pair"
  
  if(var == "flux_raw_pair")
    ret <- 
      expression(paste("Energy flux between a \nconsumer-resource pair (log J.y"^"-1"*")"))
  
  if(var == "cons_CPI_crisp")
    ret <- 
      "Crisp CPI value of consumer"
  
  if(var == "gini")
    ret <- 
      "Gini index for consumers"
  
  if(var == "CN_dist")
    ret <- 
      "Euclidean distance in CN plane\n of a consumer-resource pair"
  
  if(var == "CNP_dist")
    ret <- 
      "Euclidean distance in CNP space\n of a consumer-resource pair"
  
  if(var == "avg_CN_dist")
    ret <- 
      "Mean Euclidean distance between \nconsumers in CN plane"
  
  if(var == "avg_CNP_dist")
    ret <- 
      "Mean Euclidean distance between \nconsumers in CNP space"
  
  if(var == "totalflow_Jy")
    ret <- 
      expression(paste("Total energy flow (log J.y"^"-1"*")"))
  
  if(var == "predflow_Jy")
    ret <- 
      expression(paste("Energy flow to \npredators (log J.y"^"-1"*")"))
  
  if(var == "n_removed")
    ret <- 
      "Number of species removed \nfrom site species pool"
  
  if(var == "n_species")
    ret <- 
      "Number of species present \nin site species pool"
  
  if(var == "rsquared")
    ret <- 
      "R squared \nof model"
  
  if(var == "pd_int")
    ret <- 
      "pd of intercept"
  
  if(var == "pd_rel1")
    ret <- 
      "pd of first-order \nrelationship"
  
  if(var == "pd_rel2")
    ret <- 
      "pd of second-order \nrelationship"
  
  if(var == "rope_int")
    ret <- 
      "% in ROPE of intercept"
  
  if(var == "rope_rel1")
    ret <- 
      "% in ROPE of first-\norder relationship"
  
  if(var == "rope_rel2")
    ret <- 
      "% in ROPE of second-\norder relationship"

  # Return
  return(ret)
}

# Basic line plot ---------------------------------------------------------
basic_line <- function(model, dats, x_var, y_var, scale = "none"){
  # Get effect
  model_effect <- 
    brms::conditional_effects(model,
                              method = "fitted")[[1]] 
  # Rescale effect
  if(scale == "log"){
    model_effect <- 
      model_effect %>% 
      dplyr::mutate(dplyr::across(c(estimate__, lower__, upper__),
                                  ~exp(.)))
    # dats <- 
    #   dats %>% 
    #   dplyr::mutate(dplyr::across(dplyr::all_of(y_var),
    #                               ~log(.)))
      
  } else if(scale == "sqrt"){
    model_effect <- 
      model_effect %>% 
      dplyr::mutate(dplyr::across(c(estimate__, lower__, upper__),
                                  ~(.)^2))
    # dats <- 
    #   dats %>% 
    #   dplyr::mutate(dplyr::across(dplyr::all_of(y_var),
    #                               ~sqrt(.)))
  } else if(scale == "third"){
    model_effect <- 
      model_effect %>% 
      dplyr::mutate(dplyr::across(c(estimate__, lower__, upper__),
                                  ~(.)^3))
    # dats <- 
    #   dats %>% 
    #   dplyr::mutate(dplyr::across(dplyr::all_of(y_var),
    #                               ~sqrt(.)))
  }
  
  # Make plot
  ret <- 
    ggplot(data = model_effect,
         aes(x = .data[[x_var]],
             y = estimate__)) +
    ## Add raw data points
    geom_jitter(data = dats,
               aes(x = .data[[x_var]],
                   y = .data[[y_var]]),
               colour = "lightsteelblue4",
               alpha = 0.4) +
    ## Line
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
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) 
  
  # Return plot
  return(ret)
  
  
  
  
}


# Power and artifact line plot ---------------------------------------------------------
explo_line <- function(model, dats, x_var, y_var, scale = "none"){
 
  # Make plot
  ret <- 
    ggplot(data = dats,
           aes(x = .data[[x_var]],
               y = .data[[y_var]])) +
    geom_jitter(col = "seagreen3") +
    geom_smooth(stat = "smooth",
                col = "seagreen") +
    xlab(full_label(x_var)) +
    ylab(full_label(y_var)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"))
  
  # Return plot
  return(ret)
  
  
  
  
}



# Get euclidean distance between two points -----------
get_euclidean_distance <- function(a, b){
  
  # Get square root of sum of squares
  ret <- 
    sqrt(sum((a-b)^2))
  
  # Return
  return(ret)
}

# Get average distance between all points in a dataframe ------------------
average_dist <- function(data, CN = FALSE, CNP = FALSE){
  
  ## Remove all NAs
  data <- 
    data %>%
    dplyr::filter(!is.na(cons_C) & !is.na(cons_N) & !is.na(cons_P)) 
  
  # Get average distance
  ret <- 
    ## Get all unique, non-identical pair of points
    tidyr::crossing(id.x = 1:nrow(data), 
                    id.y = 1:nrow(data)) %>% 
    dplyr::filter(id.x != id.y) %>% 
    dplyr::filter(id.x < id.y) %>%
    ## For each one of these rows with unique combinations
    dplyr::rowwise()
  
  # Get average CN distance
  if(CN == TRUE){
    ret <- 
      ret %>% 
      ## Get the euclidean distance between the correspoding pair of points
      mutate(distance = get_euclidean_distance(c(data$cons_C[id.x], data$cons_N[id.x]),
                                               c(data$cons_C[id.y], data$cons_N[id.y]))) %>% 
      ## Ungroup to remove the rowwise
      dplyr::ungroup() %>% 
      ## Keep the distances, get mean and extract
      dplyr::select(distance) %>%
      dplyr::summarise_all(mean) %>% 
      dplyr::pull()
    
  }
  
  # Get average CNP distance
  if(CNP == TRUE){
    ret <- 
      ret %>% 
      ## Get the euclidean distance between the correspoding pair of points
      mutate(distance = get_euclidean_distance(c(data$cons_C[id.x], data$cons_N[id.x], data$cons_P[id.x]),
                                               c(data$cons_C[id.y], data$cons_N[id.y], data$cons_P[id.y]))) %>% 
      ## Ungroup to remove the rowwise
      dplyr::ungroup() %>% 
      ## Keep the distances, get mean and extract
      dplyr::select(distance) %>%
      dplyr::summarise_all(mean) %>% 
      dplyr::pull()
    
  }

  # Return
  return(ret)

}

# Functions to return list of species to remove or add ----------------------------
remove_add_specs <- function(data, number, mock_nspecs, action){
  
  
  # Select species in each
  selection <- 
    data %>% 
    dplyr::select(site, consumer) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(site) %>% 
    ## Either pick number given, or maximum number of species in species pool if lower
    slice_sample(n = ifelse(number < nrow(.), 
                            number, nrow(.)),
                 replace = F) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(consumer) %>% 
    dplyr::pull()
  
  
  # Remove species from data
  if(action == "remove"){
    ret <- 
      data %>% 
      dplyr::filter(consumer %notin% selection)} else
        # Select species from data
        if(action == "keep"){
          ret <- 
            data %>% 
            dplyr::filter(consumer %in% selection)}
  
  # Return
  return(ret)
}


# Do predictions from model from CNP --------------------------------------
do_predictions <- function(data, model_list){
  
  # Save data
  ret <- 
    tibble::tibble()
  
  
  
  # Do row by row
  for(index in 1:nrow(data)){
    print(index)
    # Get row
    row <- 
      data[index,]
    
    # Extract model
    thismodel <- 
      model_list %>% 
      dplyr::filter(taxon == row$taxon) %>% 
      dplyr::select(model)[[1]]
    
    # Blank prediction
    prediction <- 
      NA
  # No model
    if(is.na(row$nsites))
      # No model
    {
      
      prediction <- 
        tibble::tibble(Estimate = NA,
                       Est.Error = NA,
                       Q2.5 = NA,
                       Q97.5 = NA)
    } else
    
  # Model with random effect
  if(row$nsites > 1){
    prediction <- 
      row %>% 
      predict(object = thismod,
              allow_new_levels = TRUE, ### new random effects
              sample_new_levels = "gaussian") %>% 
      tibble::as_tibble()
    
  } else if (row$nsites == 1)
  # Model without random effect
  {
    prediction <- 
      data %>% 
      predict(object = thismod) %>% 
      tibble::as_tibble()
  }
    # Extract values and credible intervals
    ret <- 
      ret %>% 
      dplyr::bind_rows(dplyr::bind_cols(c(row %>% 
                                            dplyr::select(-model, -nsites), prediction)))
  }

  
  # Return
  return(ret)
  
  
  
  
}
  

# Average C, N, P values over a different taxonomic levels ----------------
average_stoich_content <- function(data, stoich){
  
  # List of levels
  levels <- 
    c("bwg_name", "genus", "family")
  
  # Data to return
  ret <- 
    data %>% 
    dplyr::mutate(C_mean = NA,
                  N_mean = NA,
                  P_mean = NA,
                  C_sd = NA,
                  N_sd = NA,
                  P_sd = NA,
                  level = NA)
  
  # Get all possible values for nutrients
  temp <-  
    tibble::tibble()
  # Average over values for each column
  for(i in levels){
    # Select column of interest and remove NAs
    temp <- 
      temp %>% 
      dplyr::bind_rows(
        stoich %>% 
        dplyr::select(do.call(paste, list(i)), C, N, P) %>% 
        dplyr::rename(taxon = 1) %>% 
        dplyr::filter(!is.na(taxon)) %>% 
        ## Now summarise
        dplyr::group_by(taxon) %>% 
        dplyr::summarise(dplyr::across(c(C, N, P),
                                       ~mean(.x, na.rm = T),
                                       .names = "{.col}_mean"),
                         dplyr::across(c(C, N, P),
                                       ~sd(.x, na.rm = T),
                                       .names = "{.col}_sd")) %>% 
        ## ungroup
        dplyr::ungroup() %>% 
        ## In case some Nan got their way in
        dplyr::mutate(dplyr::across(c(C_mean, N_mean, P_mean),
                                    ~ifelse(is.nan(.x), NA, .x))))
  }
  
  # Now go row by row, look where we have NAs, and add values
  for(i in 1:nrow(data)){
    row <- 
      data[i,]
    check <- 
      NA
    print(paste("row", i))
    ## As long as we do not have the check
    while(is.na(check)){
      for(j in levels){
        ## Extract level for which to find elemental content
        target <- 
          dplyr::pull(row[,j])
        ### In case is NA
        if(is.na(target)) {next}
        ### If not 
        val <- 
          temp %>% 
          dplyr::filter(taxon == target)
        ## Now check if exists
          ### If nothing found just leave
          if(nrow(val) == 0 & j == "ord") 
            {check <-  "y" 
             break} else
              ### If exists
              if(nrow(val) == 1){
                #### Add values to row
                ret[i,26:31] <- 
                  val[,2:7]
                ret[i,32] <- 
                  j
                #### First check is done
                check <- 
                  "y"
                break
                } else 
                  check <- 
                      "y"
      }
      
      
    }
    
    
    
    
    
  }
  
  # Return
  return(ret)
}

