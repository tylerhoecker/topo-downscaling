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
clim_vars <- c('def')

# A 'suffix' or naming convention common to all files to be downscaled
coarse_name <- 'terra'

# Time periods
times <- c('1961-1990') #,,paste0('2C_',1985:2015)

# A 'suffix' or naming convention common to the files be used as a 'template' for downscaling
templ_name <- '_topo_1981-2010.tif'

# Buffer distance for geospatial regression, in meters
buff_dist <- 50000

# Create an output directory
out_dir <- 'gids_output_nug'

if (!file.exists(paste0('../data/',out_dir))) {
  dir.create(paste0('../data/',out_dir))
}


# ------------------------------------------------------------------------------
# Run process over each tile
#-------------------------------------------------------------------------------
# Prepare inputs
#-------------------------------------------------------------------------------
# Set parrellelization plan
plan(multisession, workers = 5)
#tic()
list.files('../data/tiles') %>% 
  future_walk(\(tile){ # tile=list.files('../data/tiles')[[100]]
    
    if ( length(list.files('../data/gids_output_nug/', 
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
    
    as.list(times) %>% 
      walk(\(time){
        
        if (file.exists(paste0('../data/gids_output_nug/def_',time,tile))) {
          return(NULL)
        } else {
          # Make a log
          write(paste0('Started ',time,tile,' at ',Sys.time()),                                            
                file = "log_faster.txt",
                append = TRUE)
        }
        
        ds_coarse_tile <- rast(paste0('../data/tile_templates/ds_coarse_',time,tile))
        template_fine_tile <- rast(paste0('../data/tile_templates/template_fine_',time,tile))
        template_coarse_tile <- rast(paste0('../data/tile_templates/template_coarse_',time,tile))
        
        # Nugget distance, in meters - points within this distance are not used in regression
        # Flint and Flint suggest using the resolution of the layer to be downscaled (here, 4-km)
        nug_dist <- round(max(res(template_coarse_tile)))+1
        
        # Coarse extract for all variables, within buffer distance, returns list for each buffer
        coarse_dat <- exact_extract(ds_coarse_tile, fine_buffers, include_xy = T)
        
        # Template extract for all variables, within buffer distance, returns list for each buffer
        template_dat <- exact_extract(template_coarse_tile, fine_buffers)
        
        # Fine point extract for all variables - returns dataframe with column for each variable
        fine_dat <- terra::extract(template_fine_tile, fine_centroids, xy = T, threads = T)

        # Clean up memory
        rm(ds_coarse_tile, template_fine_tile, template_coarse_tile)
        gc()
        
        # Preallocate dataframe
        result_df <- data.frame('pt_ID' = numeric(0), 'var' = character(0), 'gids' = numeric(0))
        
        # --------------------------------------------------------------------------
        # Iterate over each point of `tile`
        #---------------------------------------------------------------------------
        focal_result_df <- imap_dfr(coarse_dat, \(pt_dat,i){
          
          # Fine-grid information from fine-grid focal location
          X <- st_coordinates(fine_centroids[i,])[,'X']
          Y <- st_coordinates(fine_centroids[i,])[,'Y']
          Cs <- fine_dat[i,]
          
          # Distance calculations for each point in this tile
          # Perhaps this could be vectorized somehow, but I couldn't figure it out
          x_dist <- X - pt_dat[['x']]
          y_dist <- Y - pt_dat[['y']]
          pt_dat[['di']] <- sqrt(x_dist^2 + y_dist^2)
          
          # Map (apply), returning dataframe, over each climate variable 
          #focal_result_df <- data.frame('pt_ID' = numeric(0), 'var' = character(0), 'gids' = numeric(0))
          
          as.list(clim_vars) |> # REPLACE HERE ALL THE YEARS... 
            map_dfr(\(clim_var){
              
              # Fine-grid information from fine-grid focal location
              C <- Cs[,clim_var]
              
              # Build dataframe with coarse-grid information from coarse grid points within buffer distance
              # Mostly this is just re-naming things so they match the Flint and Flint nomenclature
              model_df <- data.frame(
                'Zi' = pt_dat[,clim_var],
                'Xi' = pt_dat[,'x'],
                'Yi' = pt_dat[,'y'],
                'Ci' = template_dat[[i]][,clim_var],
                # Distances between fine-grid focal location and coarse centroids
                'di' = pt_dat[,'di']) 
              
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
              sum1 <- sum( (model_df$Zi + (X-model_df$Xi)*Cx + (Y-model_df$Yi)*Cy + (C-model_df$Ci)*Cc) / model_df$di)
              sum2 <- sum(1/model_df$di)
              Z = sum1/sum2 
              
              # Remove created objects that will not be needed
              rm(C, model_df, lm_mod, sum1, sum2)
              
              # Return as dataframe (foreach works like a function)
              return(data.frame('pt_ID' = i,
                                'var' = clim_var,
                                'gids' = Z))
        })
        })

        # Pivot the result to wide, for easy write-out as raster stack
        result_df_wide <- focal_result_df |> 
          pivot_wider(names_from = var, values_from = gids) |> 
          select(-pt_ID)
        
        # Create downscaled raster
        out_rast <- st_coordinates(fine_centroids) |> 
          cbind(result_df_wide) |> 
          rast(type = 'xyz', crs = tile_rast)
        
        writeRaster(out_rast, 
                    paste0('../data/gids_output_nug/',clim_vars,'_',time,tile),
                    overwrite = T)
        
        rm(focal_result_df, result_df_wide, out_rast)
        gc()
      })
  },.progress = T)

    




#toc()

# End parallelization
plan(sequential)


library(profvis)
profvis({
})




