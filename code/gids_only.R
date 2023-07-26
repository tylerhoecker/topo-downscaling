# ------------------------------------------------------------------------------
# Prepared 2/27/2023 by Tyler Hoecker: https://github.com/tylerhoecker
# ------------------------------------------------------------------------------

# This spatially downscales gridded climate data using three methods:
# - GIDS (Gradient and Inverse Distance-Squared) of Nalder and Wein (1998).
# as described in Flint and Flint (2012) and Rodman et al. (2020).
# - GLM using fine-scale data at centroid of coarse-scale data as predictor
# - GLM using aggregated fine-scale data aligned with coarse-scale as as predictor
# - All methods imploy a nugget effect.
# 
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

using('purrr','furrr','terra','geosphere','parallelly','tictoc','tools','sf')
# ------------------------------------------------------------------------------

tic()
# ------------------------------------------------------------------------------
# User-defined inputs/parameters
# ------------------------------------------------------------------------------
# Climate clim_variable - one at a time for now...
clim_var <- 'def'

# Layer names
# These are both historical datasets. 1st step is to downscale historical TerraClimate
template_suffix <- '_topo_1981-2010.tif' 
coarse_suffix <- '_terra_1981-2010.tif'

# Buffer distance for geospatial regression, in meters
buff_dist <- 50000
# Nugget distance, in meters - points within this distance are not used in regression
# Flint and Flint suggest using the resolution of the layer to be downscaled (here, 4-km)
nug_dist <- 4000

# Make tiles
# Build information about fine-scale coordinates for this tile
# Not necceasry to build tile for coarse raster
template_general <- rast(file.path('data','cogs',paste0(clim_var,template_suffix)))

# For now, crop to OR and WA
or_wa <- sf::read_sf('../../Work/GIS/cb_2018_us_state_20m/or_wa_borders.shp')
template_general <- crop(template_general, or_wa)

makeTiles(template_general, ceiling(dim(template_general)/4)[1:2], 
          filename = paste0('data/tiles/',clim_var,file_path_sans_ext(template_suffix),'_.tif'), 
          extend = TRUE, na.rm = TRUE, overwrite = TRUE)

rm(template_general)

#----------------------
# Parellelization happens over raster tiles
#----------------------
# Set parrellelization plan
plan(multisession, workers = 3)

# Run process over each tile
#tile = list.files('data/tiles')[[1]]
list.files('data/tiles') |>
  future_walk(\(tile){
    
    print(paste0('Processing tile ',tile))
    
    # Load the fine template for this tile
    template_fine <- rast(file.path('data','tiles',tile))
    
    # Load the complete
    ds_coarse <- rast(file.path('data','cogs',paste0(clim_var,coarse_suffix)))
    
    # DECISION HERE ABOUT INCLUDING CANADA
    # Crop coarse raster to focal area boundary + buffer distance to keep all
    # points that would be used
    template_border <- as.polygons(template_fine > -Inf)
    template_buffer <- buffer(template_border, width = buff_dist)
    
    ds_coarse <- crop(ds_coarse, template_buffer) 
    
    
    # Build information about fine-scale coordinates for this tile
    fine_pts <- as.points(template_fine)
    fine_coords <- crds(fine_pts)
    coarse_pts <- as.points(ds_coarse)
    coarse_coords <- crds(coarse_pts)
    
    # Buffers for all fine points, as a spatpolygon
    fine_buffers <- buffer(fine_pts, width = buff_dist)
     #names(fine_buffers) <- clim_var
    
    # Record CRS information
    out_crs <- crs(template_fine)
    
    rm(fine_pts)
    gc()
    
    # ------------------------------------------------------------------------------
    # Create a coarse (4 km) version of the "template", which will align with future data,
    # to regress against the coarse data
    
    # In the case of downscaling historical data, this means the regression will be 
    # essentially just be a bias-correction (aggregated topo historical ~ coarse terra historical), 
    # but not when downscaling projected data (aggregated downscaled terra ~ coarse terra future)
    
    # Resample, as a form of aggregation, the fine historical data to the desired coarse grid
    template_coarse <- resample(template_fine, ds_coarse, 'bilinear')
    # Maybe this saves memory? At least preserves sig figs
    template_coarse <- round(template_coarse)
    # Not sure why this is necessary... but it is
    crs(template_coarse) <- crs(ds_coarse)
    
    # Eventually, nix this
    # Range for plotting
    # if(clim_var == 'def'){plot_range = c(100,2500)} else{plot_range = c(100,2500)}
    # par(mfrow = c(3,1))
    # plot(template_fine, range = plot_range, main = 'Fine template (fine-fine)')
    # plot(ds_coarse, range = plot_range, main = 'Coarse to downscale (coarse-coarse)')
    # plot(template_coarse, range = plot_range, main = 'Aggregated fine template (coarse-fine)')
    
    # The input data that the function will assume are present in global env
    coarse_dat <- rast(list('ds_layer' = ds_coarse,
                            'template_c' = template_coarse))
    
    # Extract coarse data from within all buffers at once
    coarse_df <- terra::extract(coarse_dat, fine_buffers, cells = T, xy = T)
    
    # For Park's method, extract fine data at centroids of coarse data
    #coarse_fine_df <- terra::extract(template_fine, coarse_coords)
    
    # Turn this into a list to we can map over it (whole df too large otherwise)
    # This is slow, can't figure out how to return list from extract
    coarse_list <- split(coarse_df, ~ ID)
    #names(coarse_list) <- as.data.frame(template_fine)[,clim_var]
    
    # Before proceeding, remove all unnecessary stuff from global env
    rm(ds_coarse, template_coarse, coarse_dat, coarse_df, template_fine)
    gc()
    
    # ------------------------------------------------------------------------------
    
    downscaled_df <- coarse_list |> 
      map_dfr(\(focal_df){
        
        #-----------------------------------------------------------------------
        # Do distance calculations once, for all climate variables
        #-----------------------------------------------------------------------
        pt_id <- focal_df[1,'ID']
        
        print(paste('Running point', pt_id, 'of', length(fine_buffers), 'points'))
        
        # Measure distance between focal point and points in buffer
        dists <- c(distm(fine_coords[pt_id,],focal_df[,c('x','y')]))
        
        # Create nugget index
        buff_idx_50 <- dists > nug_dist & dists < buff_dist
        
        #-----------------------------------------------------------------------
        # From here, iterate through variables
        #-----------------------------------------------------------------------
        
        list('def','aet','tmin','tmax')
        
        # Fine-grid information from fine-grid focal location
        X <- fine_coords[pt_id,'x']
        Y <- fine_coords[pt_id,'y']
        C <- as.numeric(values(fine_buffers[pt_id]))
        
        # Build dataframe with coarse-grid information from coarse grid points within buffer distance
        model_df <- data.frame(
          'Zi' = focal_df[,'ds_layer'],
          'Xi' = focal_df[,'x'],
          'Yi' = focal_df[,'y'],
          'Ci' = focal_df[,'template_c'],
          # Distances between fine-grid focal location and coarse centroids
          # Inexing (filtering) rows with distances less than nugget
          'di' = dists)[buff_idx_50,] 
        
        # Remove NAs - happens at edges
        model_df <- model_df[!is.na(model_df[,'Ci']),]
        
        # Fit linear model
        x_lm <- model.matrix(Zi ~ Xi+Yi+Ci, data = model_df) 
        y_lm <- model_df[,'Zi']
        lm_mod <- .lm.fit(x_lm, y_lm)
        
        # Extract coefficients of linear model
        Cx <- lm_mod$coefficients[2]
        Cy <- lm_mod$coefficients[3]
        Cc <- lm_mod$coefficients[4]
        
        # Inverse distance weighting
        # This forumula is provided on page 5 of Flint and Flint 2021
        sum1 <- sum((model_df$Zi+(X-model_df$Xi)*Cx+(Y-model_df$Yi)*Cy+(C-model_df$Ci)*Cc) 
                    /model_df$di^2)
        sum2 <- sum(1/model_df$di^2)
        Z = sum1/sum2 
        
        # Return as dataframe
        return(data.frame('gids' = round(Z)))
      })
    
    fine_coords |> 
      cbind(downscaled_df[,'gids']) |> 
      rast(type = 'xyz', crs = out_crs) |> 
      writeRaster(paste0('data/gids_output/',names(coarse_pts),'_',tile), overwrite = T)
    
    })


toc()
plan(sequential)

# DELETE TILES FOR THIS clim_varIABLE!

# Testing info
complete <- rast(file.path(cogs_path,paste0('def',template_suffix)))
tile <- crop(complete, bounds)
(242.15   *(ncell(complete)/ncell(tile))) / 60 / 60

# Approximately 48 compute hours per clim_variable for entire western US  


