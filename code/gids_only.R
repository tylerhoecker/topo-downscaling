# ------------------------------------------------------------------------------
# Prepared 7/27/2023 by Tyler Hoecker: https://github.com/tylerhoecker
# ------------------------------------------------------------------------------

# # Set parrellelization plan
# plan(multisession, workers = 3)

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

using('purrr','tictoc','sf','exactextractr','terra','tidyverse','parallel','doParallel') 
# ------------------------------------------------------------------------------

# Layer names
# These are both historical datasets. 1st step is to downscale historical TerraClimate
template_suffix <- '_topo_1981-2010.tif' 
coarse_suffix <- '_terra_1981-2010.tif'

# Inputs for entire western US, do once
# Load a complete coarse template 
coarse_files <- paste0(c('def','aet','tmin','tmax'), coarse_suffix) #,'aet','tmin','tmax'
ds_coarse <- rast(file.path('data','cogs', coarse_files))
names(ds_coarse) <- c('def','aet','tmin','tmax') #,

# Create a coarse (4 km) version of the "template", which will align with future data,
# to regress against the coarse data
fine_files <- paste0(c('def','aet','tmin','tmax'), template_suffix) #,'aet','tmin','tmax'
template_fine <- rast(file.path('data','cogs', fine_files))
names(template_fine) <- c('def','aet','tmin','tmax') #,

# Resample, as a form of aggregation, the fine historical data to the desired coarse grid
template_coarse <- resample(template_fine, ds_coarse, 'bilinear')
# Not sure why this is necessary... but it is
crs(template_coarse) <- crs(ds_coarse)
names(template_coarse) <- paste0(c('def','aet','tmin','tmax'))

# Buffer distance for geospatial regression, in meters
buff_dist <- 50000
# Nugget distance, in meters - points within this distance are not used in regression
# Flint and Flint suggest using the resolution of the layer to be downscaled (here, 4-km)
nug_dist <- round(max(res(template_coarse)))+1


# Run process over each tile
# ------------------------------------------------------------------------------
list.files('data/tiles')|>
  walk(\(tile){
    print(paste0('Processing tile ',tile))
    
    # Load the fine template for this tile
    tile_rast <- rast(file.path('data','tiles',tile))
    
    fine_centroids <- st_as_sf(as.points(tile_rast))
    # Buffers for all fine points, as a spatpolygon
    # This object will increase in size as the size of the tile increases...
    fine_buffers <- st_buffer(fine_centroids, dist = buff_dist)
    
    # Coarse extract for all variables, within buffer distance, returns list for each buffer
    coarse_dat <- exact_extract(ds_coarse, fine_buffers, include_xy = T, progress = T)
    
    # Template extract for all variables, within buffer distance, returns list for each buffer
    template_dat <- exact_extract(template_coarse, fine_buffers, include_xy = T, progress = T)
 
    # Preallocate dataframe
    result_df <- data.frame('pt_ID' = numeric(0), 'var' = character(0), 'gids' = numeric(0))
    # --------------------------------------------------------------------------
    
    # Iterate over each point of `tile`
    #---------------------------------------------------------------------------
    # Best approach appears to be a loop...! Through each element of coarse and template lists, 
    # which are dataframes, and each row of focal point centroids, which is a dataframe
    
    # Set up parallelization
    parallel::detectCores()
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(
      n.cores, 
      type = "PSOCK"
    )
    #register it to be used by %dopar%
    doParallel::registerDoParallel(cl = my.cluster)
    
    # Run loop in parallel
    out <- foreach(
      i = 1:length(coarse_dat),
      .packages = c('sf','terra','purrr')
      ) %dopar% {
        
        # Distance calculations for each point in this tile
        # Perhaps this could be vectorized somehow, but I couldn't figure it out
        x_dist <- st_coordinates(fine_centroids)[i,'X'] - coarse_dat[[i]][['x']]
        y_dist <- st_coordinates(fine_centroids)[i,'Y'] - coarse_dat[[i]][['y']]
        dists <- sqrt( x_dist^2 + y_dist^2 )
        coarse_dat[[i]][['di']] <- dists
        
        # Map (apply), returning dataframe, over each climate variable 
        focal_result_df <- list('def','aet','tmin','tmax') |> # REPLACE HERE ALL THE YEARS... 
          map_dfr(\(clim_var){
            
            # Fine-grid information from fine-grid focal location
            Xp <- st_coordinates(fine_centroids[i,])[,'X']
            Yp <- st_coordinates(fine_centroids[i,])[,'Y']
            Cp <- fine_centroids[i,][['mean']]
            
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
            x_lm <- model.matrix(Zi~Xi+Yi+Ci, data = model_df) 
            y_lm <- model_df[,'Zi']
            lm_mod <- .lm.fit(x_lm, y_lm)
            
            # Save regression coefficients for GIDS formula
            Cx <- lm_mod$coefficients[2]
            Cy <- lm_mod$coefficients[3]
            Cc <- lm_mod$coefficients[4]
            
            # Inverse distance squared weighting (GIDS formula)
            # This forumula is provided on page 5 of Flint and Flint 2012
            sum1 <- sum(( model_df$Zi + (Xp-model_df$Xi)*Cx + (Yp-model_df$Yi)*Cy + (Cp-model_df$Ci)*Cc ) / model_df$di^2)
            sum2 <- sum(1/model_df$di^2)
            Z = sum1/sum2 
            
            # Return as dataframe (foreach works like a function)
            return(data.frame('pt_ID' = i,
                              'var' = clim_var,
                              'gids' = Z))
          })
      }
    
    # End parallelization
    parallel::stopCluster(cl = my.cluster)
    
    # Pivot the result to wide, for easy write-out as raster stack
    result_df_wide <- out |> 
      bind_rows() |> 
      pivot_wider(names_from = var, values_from = gids) |> 
      select(-pt_ID)
    
    # Write out each layer of stack as a tiff
    st_coordinates(fine_centroids) |> 
      cbind(result_df_wide) |> 
      rast(type = 'xyz', crs = tile_rast) |> 
      writeRaster(paste0('data/gids_output/','ds-270-m_',c('def'),tile), 
                  overwrite = T)
    toc()
  })




