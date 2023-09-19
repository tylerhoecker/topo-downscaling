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

using('purrr','sf','exactextractr','terra','tidyverse','furrr') 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set some parameters, to make it easier to apply to different projects
# ------------------------------------------------------------------------------
# The climate variable short names as they appear in filenames
clim_vars <- c('aet','def','tmax','tmin')

# A 'suffix' or naming convention common to all files to be downscaled
coarse_name <- 'terra'

# Time periods
times <- c('1961-1990','2C_1985-2015') #,,paste0('2C_',1985:2015)

# A 'suffix' or naming convention common to the files be used as a 'template' for downscaling
templ_name <- '_topo_1981-2010.tif'

# Buffer distance for geospatial regression, in meters
buff_dist <- 50000

# Create an output directory
out_dir <- 'gids_output'

if (!file.exists(paste0('../data/',out_dir))) {
  dir.create(paste0('../data/',out_dir))
}

# ------------------------------------------------------------------------------
# Run process over each tile
#-------------------------------------------------------------------------------
# Prepare inputs
#-------------------------------------------------------------------------------
# Set parrellelization plan
plan(multisession, workers = 20)
#tic()
list.files('../data/tiles') %>% 
  future_walk(\(tile){ # tile=list.files('../data/tiles')[[1]]
    
    if ( length(list.files('../data/gids_output/', 
                           pattern = paste0('.*',tile)))
         == length(clim_vars)*length(times)) {
      return(NULL)
    }
    
    # Load the fine template for this tile
    tile_rast <- rast(paste0('../data/tiles/',tile))
    
    # Collect the centroids of each grid cell in the tile 
    fine_centroids <- st_as_sf(as.points(tile_rast))

    # Buffers for all fine points, as a spatpolygon
    fine_buffers <- st_buffer(fine_centroids, dist = buff_dist)
    
    # Fine template is the same for all time periods, load here
    template_fine <- rast(paste0('../data/cogs/',clim_vars,templ_name))
    names(template_fine) <- clim_vars
    
    # --
    # Resample, as a form of aggregation, the fine data to coarse grid, so they can be regressed
    # Need to load one of the coarse files for this step
    eg_coarse <- rast(paste0('../data/cogs/',clim_vars[1],'_',coarse_name,'_',times[1],'.tif'))
    template_coarse <- resample(template_fine, eg_coarse, 'bilinear')
    
    rm(eg_coarse)
    gc()
    
    # Nugget distance, in meters - points within this distance are not used in regression
    # Flint and Flint suggest using the resolution of the layer to be downscaled (here, 4-km)
    nug_dist <- round(max(res(template_coarse)))+1
    
    as.list(times) %>% 
      walk(\(time){
        
        if (file.exists(paste0('../data/gids_output/def_',time,tile))) {
          return(NULL)
        } else {
          # Make a log
          write(paste0('Started ',time,tile,' at ',Sys.time()),                                            
                file = "log_faster.txt",
                append = TRUE)
        }
        
        # Stack of coarse files (all variables) for this time --
        coarse_files <- paste0('../data/cogs/',clim_vars,'_',coarse_name,'_',time,'.tif')
        ds_coarse <- rast(coarse_files)
        names(ds_coarse) <- clim_vars
        
        # Coarse extract for all variables, within buffer distance, returns list for each buffer
        coarse_dat <- exact_extract(ds_coarse, fine_buffers, include_xy = T)
        
        # Template extract for all variables, within buffer distance, returns list for each buffer
        template_dat <- exact_extract(template_coarse, fine_buffers, include_xy = T)
        
        # Fine point extract for all variables - returns dataframe with column for each variable
        fine_dat <- terra::extract(template_fine, fine_centroids, xy = T)
        
        # Clean up memory
        rm(fine_buffers, coarse_files, ds_coarse, template_coarse, template_fine)
        gc()
        
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
          dists <- sqrt(x_dist^2 + y_dist^2)
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
              sum1 <- sum( (model_df$Zi + (X-model_df$Xi)*Cx + (Y-model_df$Yi)*Cy + (C-model_df$Ci)*Cc) / model_df$di^2)
              sum2 <- sum(1/model_df$di^2)
              Z = sum1/sum2 
              
              # Remove created objects that will not be needed
              rm(fine_centroids, coarse_dat, template_dat, fine_dat, X, Y, C, model_df, lm_mod, sum1, sum2)
              
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
        st_coordinates(fine_centroids) |> 
          cbind(result_df_wide) |> 
          rast(type = 'xyz', crs = tile_rast) %>% 
          writeRaster(paste0('../data/gids_output/',clim_vars,'_',time,tile),
                      overwrite = T)
        
        rm(focal_result_df, result_df, result_df_wide)
        gc()
      })
  },.progress = T)


#toc()

# End parallelization
plan(sequential)



library(profvis)
profvis({
})




