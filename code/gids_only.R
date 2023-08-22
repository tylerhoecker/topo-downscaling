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

using('purrr','furrr','terra','geosphere','parallelly','tictoc','tools','sf')
# ------------------------------------------------------------------------------

# tic()
# ------------------------------------------------------------------------------
# User-defined inputs/parameters
# ------------------------------------------------------------------------------
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
# Not necessary to build tile for coarse raster
# Climate variable for this step is arbitary - spatial info identical for all
template_general <- rast(file.path('data','cogs',paste0('def',template_suffix)))

# For now, crop to OR and WA
or_wa <- sf::read_sf('../../Work/GIS/cb_2018_us_state_20m/or_wa_borders.shp') |> 
  st_transform("EPSG:3741")
template_general <- crop(template_general, or_wa)
rm(or_wa)
makeTiles(template_general, ceiling(dim(template_general)/5)[1:2], 
          filename = paste0('data/tiles/',file_path_sans_ext(template_suffix),'_.tif'), 
          extend = TRUE, na.rm = TRUE, overwrite = TRUE)
rm(template_general)

# Run process over each tile
# ------------------------------------------------------------------------------
# RUN THIS: 
tile = list.files('data/tiles')[[2]]

#list.files('data/tiles') |>
  #walk(\(tile){
    
    print(paste0('Processing tile ',tile))
    
    # Load the fine template for this tile
    master_fine <- rast(file.path('data','tiles',tile))
    
    # Load a complete coarse template 
    coarse_files <- paste0(c('def','aet','tmin','tmax'), coarse_suffix) #,'aet','tmin','tmax'
    ds_coarse_stack <- rast(file.path('data','cogs', coarse_files))
    names(ds_coarse_stack) <- c('def','aet','tmin','tmax') #,
    # Don't crop the coarse stack here. It's coarse, so not very heavy, and 
    # info will be needed from beyond the bounds of the tile (+ 50 km)
    # Could buffer the tile and then crop, seems unneccesary. 
    
    
    # Buffers for all fine points, as a spatpolygon
    # This object will increase in size as the size of the tile increases...
    fine_buffers <- buffer(as.points(master_fine), width = buff_dist)
    
    # Preallocate dataframe
    result_df <- data.frame('var' = character(0), 'gids' = numeric(0))
    
    # --------------------------------------------------------------------------
    # focal_buff = fine_buffers[1]
    tic()
    result_list <- split(fine_buffers[1], 'mean') |> 
      map_dfr(\(focal_buff){
        
        print(paste('Running point', crds(centroids(focal_buff))[1],crds(centroids(focal_buff))[2]))
        
        # clim_var = 'tmax'
        downscaled_df <- list('def','aet','tmin','tmax') |> # REPLACE HERE ALL THE YEARS... 
          map_dfr(\(clim_var){
            
            print(paste('Running',clim_var))
            # Create a coarse (4 km) version of the "template", which will align with future data,
            # to regress against the coarse data
            
            template_fine <- rast(file.path('data','cogs',paste0(clim_var,template_suffix)))
            #template_fine <- crop(template_fine, master_fine)
            
            # Resample, as a form of aggregation, the fine historical data to the desired coarse grid
            template_coarse <- resample(template_fine, ds_coarse_stack, 'bilinear')
            #template_coarse <- crop(template_coarse, master_fine)
            # Maybe this saves memory? At least preserves sig figs
            #template_coarse <- round(template_coarse)
            # Not sure why this is necessary... but it is
            crs(template_coarse) <- crs(ds_coarse_stack)
            names(template_coarse) <- 'template_c'
            
            # Extract...
            coarse_dat <- terra::extract(ds_coarse_stack[[clim_var]], focal_buff, xy = T)
            template_dat <- terra::extract(template_coarse, focal_buff, xy = T)
            # Combine, for convenience
            focal_df <- cbind(coarse_dat, 'template_c' = template_dat[['template_c']])
            
            # Create spatial vector of coordinates of cells within buffer
            focal_vect <- vect(focal_df, geom=c("x", "y"), crs = crs(focal_buff))
            
            # Measure distance between focal point and points/cels in buffer
            dists <- c(distance(centroids(focal_buff), focal_vect))
            
            # Fine-grid information from fine-grid focal location
            X <- crds(centroids(focal_buff))[,'x']
            Y <- crds(centroids(focal_buff))[,'y']
            C <- terra::extract(template_fine, crds(centroids(focal_buff)))[['mean']]
            
            # Build dataframe with coarse-grid information from coarse grid points within buffer distance
            # Mostly this is just re-naming things so they match the Flint and Flint nomenclature
            model_df <- data.frame(
              'Zi' = focal_df[,clim_var],
              'Xi' = focal_df[,'x'],
              'Yi' = focal_df[,'y'],
              'Ci' = focal_df[,'template_c'],
              # Distances between fine-grid focal location and coarse centroids
              'di' = dists) 
            
            # Remove NAs - happens at edges
            model_df <- na.omit(model_df)
            
            # Nugget effect - remove points within distance of coarse cell resolution (4000)
            # Eliminates possible halo effect when fine and coarse cell centroids are near
            model_df <- model_df[model_df[,'di'] > 4000,] 
            
            # Fit linear model
            #lm_mod <-lm (Zi ~ Xi + Yi + Ci, data = model_df) 
            x_lm <- model.matrix(Zi ~ Xi + Yi + Ci, data = model_df)
            y_lm <- model_df[,'Zi']
            lm_mod <- .lm.fit(x_lm, y_lm)

            # Extract coefficients of linear model
            Cx <- lm_mod$coefficients[2]
            Cy <- lm_mod$coefficients[3]
            Cc <- lm_mod$coefficients[4]
            
            # Inverse distance weighting
            # This forumula is provided on page 5 of Flint and Flint 2012
            sum1 <- sum(( model_df$Zi + (X-model_df$Xi)*Cx + (Y-model_df$Yi)*Cy + (C-model_df$Ci)*Cc ) / model_df$di^2)
            sum2 <- sum(1/model_df$di^2)
            Z = sum1/sum2 
            
            # Return as dataframe
            return(data.frame('var' = clim_var,
                              'gids' = Z))
          })
        # End map_dfr(\(clim_var)  
        result_df <- rbind(result_df, downscaled_df)
        return(result_df)
      }, .id = 'cell_ID')
    toc()
    # End walk(\(focal_buff)

    
    
          
            
              
              crds(centroids(fine_buffers)) |> 
                cbind(downscaled_df[,'gids']) |> 
                rast(type = 'xyz', crs = out_crs) |> 
                writeRaster(paste0('data/gids_output/',names(coarse_pts),'_',tile), overwrite = T)
              
              crds(centroids(fine_buffers)) |> 
                cbind(downscaled_df[,'glm']) |> 
                rast(type = 'xyz', crs = out_crs) |> 
                writeRaster(paste0('data/glm_output/',names(coarse_pts),'_',tile), overwrite = T)




toc()
plan(sequential)

# DELETE TILES FOR THIS clim_varIABLE!

# Testing info
complete <- rast(file.path(cogs_path,paste0('def',template_suffix)))
tile <- crop(complete, bounds)
(242.15   *(ncell(complete)/ncell(tile))) / 60 / 60

# Approximately 48 compute hours per clim_variable for entire western US  


