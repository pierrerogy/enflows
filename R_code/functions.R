# Functions
# Author: Pierre Rogy

# Libraries
library(tidyverse)
library(curl)
library(here)
library(raster)
library(sp)
library(cheddar)

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
        dplyr::pull() %>% 
        ## Add basal species
        vctrs::vec_c(basal)
      
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
    mat %>% 
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
  