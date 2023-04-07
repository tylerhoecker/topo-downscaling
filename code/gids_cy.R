# ------------------------------------------------------------------------------
# Prepared 2/27/2023 by Tyler Hoecker: https://github.com/tylerhoecker
#
# Adapted from Rodman, K C., et al. 2020. 
# A Changing Climate is Snuffing Out Post-Fire Recovery in Montane Forests. 
# Global Ecology and Biogeography.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
## This code performs a multi-step procedure that spatially downscales gridded climate
# data using GIDS (Gradient and Inverse Distance-Squared) of Nalder and Wein (1998)
# as described in Flint and Flint (2012)
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

using('purrr','furrr','terra','geosphere','parallelly','tictoc')
# ------------------------------------------------------------------------------

tic()
# ------------------------------------------------------------------------------
# User-defined inputs/parameters
# ------------------------------------------------------------------------------
# For now, use local file paths...
cogs_path <- file.path('data','cogs') # Use correct paths for CyVerse: cogs_path <- file.path('analyses','gids','cogs')

# Layer names
# These are both historical datasets. 1st step is to downscale historical TerraClimate
template_suffix <- '_topo_1981_2010_cog.tif' 
ds_suffix <- '_hist_1981_2010_cog.tif'

# Optional bounding box - nix eventually as this will iterate on tiles
bounds <- ext(-112.5, -111.5, 36, 36.5) # or NA 

# Buffer distance for geospatial regression, in meters
buff_dist <- 15000
# Nugget distance, in meters - points within this distance are not used in regression
# Flint and Flint suggest using the resolution of the layer to be downscaled (here, 4-km)
nug_dist <- 4000

#----------------------
# Ultimately, parellelization happens over raster tiles
#----------------------
# Set parrellelization plan
#plan(multisession, workers = availableCores()-2)

# Build information about fine-scale coordinates for this tile
template_general <- rast(file.path(cogs_path,paste0('def',template_suffix)))
coarse_general <- rast(file.path(cogs_path,paste0('def',ds_suffix)))
  
if(!is.na(bounds)){
  template_general <- crop(template_general, bounds)
  coarse_general <- crop(coarse_general, bounds)
}

# As points (centroids)
fine_pts <- as.points(template_general)
fine_coords <- crds(fine_pts)
coarse_pts <- as.points(coarse_general)
coarse_coords <- crds(coarse_pts)

# Buffers for all fine points, as a spatpolygon
fine_buffers <- buffer(fine_pts, width = buff_dist)
#names(fine_buffers) <- var

# Save crs
out_crs <- crs(template_general)


# ------------------------------------------------------------------------------
# var = 'def'
'def' |>
  walk(function(var){
    
    # Fine-scale template 
    template_fine <- rast(file.path(cogs_path,paste0(var,template_suffix)))
    names(template_fine) <- var
    
    # Raster that will be downscaled
    ds_coarse <- rast(file.path(cogs_path,paste0(var,ds_suffix)))

    # Crop rasters to focal area boundary, if provided
    if(!is.na(bounds)){
      template_fine <- crop(template_fine, bounds)
      ds_coarse <- crop(ds_coarse, bounds) 
      # Terra climate also needs to be masked because it extends into Mexico/Canada and that will 
      # change the points used in the regression, but extents don't match... use CONUS polygon?
    }
    
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
    crs(template_coarse) <- ds_coarse
    
    # Eventually, nix this
    # Range for plotting
    if(var == 'def'){plot_range = c(100,2500)} else{plot_range = c(100,2500)}
    par(mfrow = c(3,1))
    plot(template_fine, range = plot_range, main = 'Fine template (fine-fine)')
    plot(ds_coarse, range = plot_range, main = 'Coarse to downscale (coarse-coarse)')
    plot(template_coarse, range = plot_range, main = 'Aggregated fine template (coarse-fine)')
    
    # The input data that the function will assume are present in global env
    coarse_dat <- rast(list('ds_layer' = ds_coarse,
                            'template_c' = template_coarse))
    
    # Extract coarse data from within all buffers at once
    coarse_df <- terra::extract(coarse_dat, fine_buffers, cells = T, xy = T)
    
    # For Park's method, extract fine data at centroids of coarse data
    coarse_fine_df <- terra::extract(template_fine, coarse_coords)
    
    # Turn this into a list to we can map over it (whole df too large otherwise)
    # This is slow, can't figure out how to return list from extract
    coarse_list <- split(coarse_df, ~ ID)
    names(coarse_list) <- as.data.frame(template_fine)[,var]
    
   # Before proceeding, remove all unnecessary stuff from global env
    rm(ds_coarse, template_general, template_coarse, coarse_dat, coarse_general, 
       coarse_df, template_fine, fine_pts)
    gc()
    
    # ------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------
    downscaled_df <- coarse_list |> 
      # imap syntax: focal_df = element i of list of dataframes, focal_idx = index (name of element), which is the fine value
      imap_dfr(\(focal_df,focal_idx){
        
        pt_id <- focal_df[1,'ID']
        
        print(paste('Running point', pt_id, 'of', length(fine_buffers), 'points'))
        
        # Measure distance between focal point and points in buffer
        dists <- c(distm(fine_coords[pt_id,],focal_df[,c('x','y')]))
        
        # Create nugget index
        nug_idx <- dists > nug_dist
        
        # Fine-grid information from fine-grid focal location
        X <- fine_coords[pt_id,'x']
        Y <- fine_coords[pt_id,'y']
        E <- as.numeric(focal_idx)
        
        # Build dataframe with coarse-grid information from coarse grid points within buffer distance
        model_df <- data.frame(
          'Zi' = focal_df[,'ds_layer'],
          'Xi' = focal_df[,'x'],
          'Yi' = focal_df[,'y'],
          'Ei' = focal_df[,'template_c'],
          'Pi' = coarse_fine_df[focal_df[,'cell'],],
          # Distances between fine-grid focal location and coarse centroids
          # Inexing (filtering) rows with distances less than nugget
          'di' = dists)[nug_idx,]
        
        # Fit linear model
        x_lm <- model.matrix(Zi ~ Xi+Yi+Ei, data = model_df) 
        y_lm <- model_df[,'Zi']
        lm_mod <- .lm.fit(x_lm, y_lm)
        
        # Extract coefficients of linear model
        Cx <- lm_mod$coefficients[2]
        Cy <- lm_mod$coefficients[3]
        Ce <- lm_mod$coefficients[4]
        
        # Inverse distance weighting
        # This forumula is provided on page 5 of Flint and Flint 2021
        sum1 <- sum((model_df$Zi+(X-model_df$Xi)*Cx+(Y-model_df$Yi)*Cy+(E-model_df$Ei)*Ce) 
                    /model_df$di^2)
        sum2 <- sum(1/model_df$di^2)
        Z = sum1/sum2 # 1016.312

        # ----------------------------------------------------------------------
        # Parks version
        # ----------------------------------------------------------------------
        # Apparently predict doesn't fork on fit.glm, so using the normal way instead, may be slower
        glm_mod_agg <- glm(Zi ~ Ei, data = model_df, family='quasipoisson')
        glm_pred_agg <- predict(glm_mod_agg, newdata = data.frame('Ei' = E), type='response')
        
        glm_mod_fine <- glm(Zi ~ Pi, data = model_df, family='quasipoisson')
        glm_pred_fine <- predict(glm_mod_fine, newdata = data.frame('Pi' = E), type='response')
        
        
        # Return as dataframe
        return(data.frame('gids' = round(Z),
                          'glm_agg' = round(glm_pred_agg),
                          'glm_fine' = round(glm_pred_fine)))
      })
  
    # Can shorten this eventually... but whatever   
    # Write raster for each method
    gids_df <- cbind(fine_coords, downscaled_df[,'gids'])
    gids_rast <- rast(gids_df, type = 'xyz', crs = out_crs)
    writeRaster(gids_rast, paste0('output/def/grand_gids.tiff'), overwrite = T)
    
    glm_agg_df <- cbind(fine_coords, downscaled_df[,'glm_agg'])
    glm_agg_rast <- rast(glm_agg_df, type = 'xyz', crs = out_crs)
    writeRaster(glm_agg_rast, paste0('output/def/grand_glm_agg.tiff'), overwrite = T)
    
    glm_fine_df <- cbind(fine_coords, downscaled_df[,'glm_fine'])
    glm_fine_rast <- rast(glm_fine_df, type = 'xyz', crs = out_crs)
    writeRaster(glm_fine_rast, paste0('output/def/grand_glm_fine.tiff'), overwrite = T)
    
  })
toc()

359*ncell(rast(file.path(cogs_path,paste0('def',template_suffix))))/ncell(template_general) /
  60 / 60

# Approximately 80 compute hours per variable for entire western US  


