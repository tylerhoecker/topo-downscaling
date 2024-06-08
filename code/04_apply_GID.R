# ------------------------------------------------------------------------------
# Prepared by Tyler Hoecker: https://github.com/tylerhoecker
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Import  packages
# ------------------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(purrr)
library(furrr)
library(sf)
library(exactextractr)
library(terra)
# library(future.callr)
library(data.table)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Set some parameters, to make it easier to apply to different projects
# ------------------------------------------------------------------------------

# The climate variable short names as they appear in filenames
clim_vars <- c('aet','def','tmax','tmin')

# A 'suffix' or naming convention common to all files to be downscaled
coarse_name <- 'terra_2C'
coarse_name <- 'terra_hist'
# Time periods
times <- c(1985:2015) #paste0('2C_',1985:2015) #c('1961-1990','2C_1985-2015') 
times <- c(1961:2022)
# A 'suffix' or naming convention common to the files be used as a 'template' for downscaling
templ_name <- '_topo_1981-2010.tif'

# Buffer distance for geospatial regression, in meters
buff_dist <- 50000

# Create a data directory
data_root <- 'data/'

# Create an output directory
out_dir <- 'gids_output/'

if (!file.exists(paste0(data_root,out_dir))) {
  dir.create(paste0(data_root,out_dir))
}

# Directory for tile_templates (cropped rasters for each tile)
tile_templ_dir <- 'data/tile_templates/'

# Find incomplete tiles
done <- times %>% 
  map(\(time){
    list.files(paste0(data_root,out_dir), pattern = paste0(time,'_.*'))
  }) %>% 
  list_c()

to_do <- times %>% 
  map(\(time){
    paste0(time,list.files(paste0(data_root,'tiles/')))
  }) %>% 
  list_c()

not_done <- unique(sub('.*_','_',to_do[!to_do %in% done]))


# ------------------------------------------------------------------------------
# Run process over each tile
#-------------------------------------------------------------------------------
# Prepare inputs
#-------------------------------------------------------------------------------
# Set parrellelization plan
availableCores()
 plan(
   list(
     future::tweak(
       multisession,
       workers = (availableCores() - 2) %/% 6),
     future::tweak(
       multisession,
       workers = 6)
     )
   )
# plan(multisession, workers = 14)
options(mc.cores = 16)
not_done %>% 
  future_walk(\(tile){ # tile=not_done[[1]]
    options(mc.cores =6)
    if ( length(list.files(paste0(data_root,out_dir), 
                           pattern = paste0('.*',tile)))
         == length(times)) {
      return(NULL)
    }
    
    # Load the fine template for this tile
    tile_rast <- rast(paste0(data_root,'tiles/',tile))
    # Save the CRS info
    crs_out <- crs(tile_rast)
    
    # Collect the centroids of each grid cell in the tile 
    fine_centroids <- as.points(tile_rast)
    # Reproject these points to UTM, then cbind()
    
    # Buffers for all fine points, as a spatpolygon
    fine_buffers <- st_as_sf(buffer(fine_centroids, buff_dist))
    
    #load in templates
    template_fine <- rast(paste0(tile_templ_dir,'template_fine_.tif'))
    template_coarse <- rast(paste0(tile_templ_dir,'template_coarse_.tif'))
    
    as.list(times) %>% 
      walk(\(time){
        tic()
        if (file.exists(paste0(data_root,out_dir,coarse_name,"_",time,tile))) {
          return(NULL)
        } else {
          # Make a log
          write(paste0('Started ',time,tile,' at ',Sys.time()),                                            
                file = "log_faster.txt",
                append = TRUE)
        }
        
        # Skip tiles that have not been made yet
        if( file.exists(paste0(tile_templ_dir,'ds_coarse_',coarse_name,'_',time,'.tif')) ){
          ds_coarse <- rast(paste0(tile_templ_dir,'ds_coarse_',coarse_name,'_',time,'.tif'))
          
          
        } else {
          return(NULL)
        }
        
        # Nugget distance, in meters - points within this distance are not used in regression
        # Flint and Flint suggest using the resolution of the layer to be downscaled (here, 4-km)
        nug_dist <- 4000 / 111000
        
        # Coarse extract for all variables, within buffer distance, returns list for each buffer
        coarse_dat <- exact_extract(ds_coarse, fine_buffers, include_xy = T)
        
        # Template extract for all variables, within buffer distance, returns list for each buffer
        template_dat <- exact_extract(template_coarse, fine_buffers)
        
        # Fine point extract for all variables - returns dataframe with column for each variable
        fine_dat <- terra::extract(template_fine, fine_centroids, xy = T, threads = T)
        
        # Clean up memory
        rm(ds_coarse)
        #gc()
        toc()
        # --------------------------------------------------------------------------
        # Iterate over each point of `tile`
        #---------------------------------------------------------------------------
        tic()
        focal_result_df <- future_imap_dfr(coarse_dat, \(pt_dat,i){
          
          # Fine-grid information from fine-grid focal location
          X <- fine_dat[i,'x']
          Y <- fine_dat[i,'y']
          Cs <- fine_dat[i,]
          
          # Distance calculations for each point in this tile
          # Perhaps this could be vectorized somehow, but I couldn't figure it out
          x_dist <- X - pt_dat[['x']]
          y_dist <- Y - pt_dat[['y']]
          pt_dat[['di']] <- sqrt(x_dist^2 + y_dist^2)
          
          # Map (apply), returning dataframe, over each climate variable 
          #focal_result_df <- data.frame('pt_ID' = numeric(0), 'var' = character(0), 'gids' = numeric(0))
          
          as.list(clim_vars) %>%  # REPLACE HERE ALL THE YEARS... 
            map_dfr(\(clim_var){
              
              # Fine-grid information from fine-grid focal location
              C <- Cs[,clim_var]
              
              # Build dataframe with coarse-grid information from coarse grid points within buffer distance
              # Mostly this is just re-naming things so they match the Flint and Flint nomenclature
              model_df <- data.table(
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
              model_df <- model_df[di > nug_dist] 
              
              # The regression
              # lm_mod <- lm(Zi ~ Xi + Yi + Ci, data = model_df) 
              x_lm <- model.matrix(Zi~Xi+Yi+Ci, data = model_df) 
              y_lm <- as.matrix(model_df[, c('Zi')])
              lm_mod <- .lm.fit(x_lm, y_lm)
              
              # Save regression coefficients for GIDS formula
              Cx <- lm_mod$coefficients[2]
              Cy <- lm_mod$coefficients[3]
              Cc <- lm_mod$coefficients[4]
              
              # Inverse distance squared weighting (GIDS formula)
              # This formula is provided on page 5 of Flint and Flint 2012
              # Formula was revised - distance are not squared here - to reduce circular artifacts in areas with low 
              # variability in climate or near edges
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
        if ( length(unique(fine_dat[,'x'])) == 1 | length(unique(fine_dat[,'y'])) == 1 ){
         tile_rast <- rast(rep(paste0(data_root,'tiles/',tile),4))
         out_rast <- setValues(tile_rast, result_df_wide) 
         out_rast <- out_rast*(tile_rast/tile_rast)
        } else {
          out_rast <-  fine_dat[,c('x','y')] |> 
            cbind(result_df_wide) |> 
            rast(type = 'xyz', crs = crs_out)
        }
        toc()
        writeRaster(out_rast, 
                    paste0(data_root,out_dir,coarse_name,"_",time,tile),
                    overwrite = T)
        
        rm(focal_result_df, result_df_wide, out_rast)
        gc()
      })
  },.progress = T)

  plan(sequential)
