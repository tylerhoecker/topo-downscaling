# ------------------------------------------------------------------------------
# Prepared by Tyler Hoecker: https://github.com/tylerhoecker
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Import  packages
# ------------------------------------------------------------------------------
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs,require,character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

using('purrr','tictoc','sf','exactextractr','terra','tidyverse','furrr') 
# ------------------------------------------------------------------------------
setwd('C:/Users/PC/Desktop/tyler_working')

# Run process over each tile
# ------------------------------------------------------------------------------
# Set parrellelization plan
plan(multisession, workers = 10)

# The variables names
clim_vars <- c('def','aet','tmin','tmax')
# Buffer distance for geospatial regression, in meters
buff_dist <- 50000

tic('Total run time 5 tiles, 10 core')

list.files('data/tiles')[3010:3015] %>% 
  future_walk(\(tile){
    
    if(file.exists(paste0('data/gids_output/ds-220-m_','def',tile))){
      return(NULL)
    }
    
    # Message
    # print(paste0('Processing tile ',tile))
    # Make a log
    write(paste0('Started tile ',tile,'at ',Sys.time()),                                            
          file = "log_faster.txt",
          append = TRUE)
    
    # Load template and to-downscale layers that were prepared for entire
    # study area in `make_tiles` script
    # ***Create this for each time period of interest...***
    ds_coarse <- rast(paste0('data/tile_templates/ds_coarse',tile)) 
    template_coarse <- rast(paste0('data/tile_templates/template_coarse',tile))
    template_fine <- rast(paste0('data/tile_templates/template_fine',tile))
    
    # Nugget distance, in meters - points within this distance are not used in regression
    # Flint and Flint suggest using the resolution of the layer to be downscaled (here, 4-km)
    nug_dist <- round(max(res(template_coarse)))+1
    
    # Load the fine template for this tile
    tile_rast <- rast(file.path('data','tiles',tile))
    
    fine_centroids <- st_as_sf(as.points(tile_rast))
    # Buffers for all fine points, as a spatpolygon
    # This object will increase in size as the size of the tile increases...
    fine_buffers <- st_buffer(fine_centroids, dist = buff_dist)
    
    # Coarse extract for all variables, within buffer distance, returns list for each buffer
    coarse_dat <- exact_extract(ds_coarse, fine_buffers, include_xy = T)
    
    # Template extract for all variables, within buffer distance, returns list for each buffer
    template_dat <- exact_extract(template_coarse, fine_buffers, include_xy = T)
    
    # Fine point extract for all variables - returns dataframe with column for each variable
    fine_dat <- terra::extract(template_fine, fine_centroids, xy = T)
    
    # Clean up memory
    rm(ds_coarse, template_fine, template_coarse)

    # Preallocate dataframe
    result_df <- data.frame('pt_ID' = numeric(0), 'var' = character(0), 'gids' = numeric(0))
    # --------------------------------------------------------------------------
    
    # Iterate over each point of `tile`
    #---------------------------------------------------------------------------
    # Best approach appears to be a loop...! Through each element of coarse and template lists, 
    # which are dataframes, and each row of focal point centroids, which is a dataframe
    for (i in 1:length(coarse_dat)){    
      
      # Distance calculations for each point in this tile
      # Perhaps this could be vectorized somehow, but I couldn't figure it out
      x_dist <- st_coordinates(fine_centroids)[i,'X'] - coarse_dat[[i]][['x']]
      y_dist <- st_coordinates(fine_centroids)[i,'Y'] - coarse_dat[[i]][['y']]
      dists <- sqrt( x_dist^2 + y_dist^2 )
      coarse_dat[[i]][['di']] <- dists
      
      # Map (apply), returning dataframe, over each climate variable 
      focal_result_df <- data.frame('pt_ID' = numeric(0), 'var' = character(0), 'gids' = numeric(0))
      focal_result_df <- as.list(clim_vars) |> # REPLACE HERE ALL THE YEARS... 
        map_dfr(\(clim_var){
          
          # Fine-grid information from fine-grid focal location
          X <- st_coordinates(fine_centroids[i,])[,'X']
          Y <- st_coordinates(fine_centroids[i,])[,'Y']
          C <- fine_dat[i,clim_var]
          
          # Build dataframe with coarse-grid information from coarse grid points within buffer distance
          # Mostly this is just re-naming things so they match the Flint and Flint nomenclature
          model_df <- data.frame(
            'Zi' = coarse_dat[[i]][,clim_var],
            'Xi' = coarse_dat[[i]][,'x'],
            'Yi' = coarse_dat[[i]][,'y'],
            'Ci' = template_dat[[i]][,clim_var],
            # Distances between fine-grid focal location and coarse centroids
            'di' = coarse_dat[[i]][,'di']) 
          
          # Remove NAs - happens at edges
          model_df <- na.omit(model_df)
          
          # Nugget effect - remove points within distance of coarse cell resolution (4000)
          # Eliminates possible halo effect when fine and coarse cell centroids are near
          model_df <- model_df[model_df[,'di'] > nug_dist,] 
          
          # The regression
          # lm_mod <- lm(Zi ~ Xi + Yi + Ci, data = model_df) 
          x_lm <- model.matrix(~Xi+Yi+Ci, data = model_df) 
          y_lm <- model_df[,'Zi']
          lm_mod <- .lm.fit(x_lm, y_lm)
          
          # Save regression coefficients for GIDS formula
          Cx <- lm_mod$coefficients[2]
          Cy <- lm_mod$coefficients[3]
          Cc <- lm_mod$coefficients[4]
          
          # Inverse distance squared weighting (GIDS formula)
          # This formula is provided on page 5 of Flint and Flint 2012
          sum1 <- sum(( model_df$Zi + (X-model_df$Xi)*Cx + (Y-model_df$Yi)*Cy + (C-model_df$Ci)*Cc ) / model_df$di^2)
          sum2 <- sum(1/model_df$di^2)
          Z = sum1/sum2 
          
          # Return as dataframe (foreach works like a function)
          return(data.frame('pt_ID' = i,
                            'var' = clim_var,
                            'gids' = Z))
        })
      result_df <- rbind(result_df, focal_result_df)
    }
    
    # Pivot the result to wide, for easy write-out as raster stack
    result_df_wide <- result_df |> 
      pivot_wider(names_from = var, values_from = gids) |> 
      select(-pt_ID)
    
    # Create downscaled raster
    out_rast <- st_coordinates(fine_centroids) |> 
      cbind(result_df_wide) |> 
      rast(type = 'xyz', crs = tile_rast) %>% 
      writeRaster(paste0('data/gids_output/ds-220-m_',clim_vars,tile),
                  overwrite = T)
    
    rm(focal_result_df, result_df, result_df_wide, tile_rast, 
       fine_centroids, fine_buffers, coarse_dat, template_dat, fine_dat)
    gc()
  })

toc()
# End parallelization
plan(sequential)


  



