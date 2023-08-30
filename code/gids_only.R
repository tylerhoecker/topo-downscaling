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

using('purrr','furrr','tictoc','sf','exactextractr','terra','tools','tidyverse','svMisc','Rcpp','RcppArmadillo') #'terra','geosphere','parallelly','tools','stars'

# Compile C++ version of linear model 
# ------------------------------------------------------------------------------
Rcpp::sourceCpp("code/lm_rcpp.cpp")

# ------------------------------------------------------------------------------

# tic()
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
# Layer names
# These are both historical datasets. 1st step is to downscale historical TerraClimate
template_suffix <- '_topo_1981-2010.tif' 
coarse_suffix <- '_terra_1981-2010.tif'

# Make tiles -------------------------------------------------------------------
# Build information about fine-scale coordinates for this tile
# Not necessary to build tile for coarse raster
# Climate variable for this step is arbitary - spatial info identical for all
template_general <- rast(file.path('data','cogs',paste0('def',template_suffix)))

# For now, crop to OR and WA
# or_wa <- sf::read_sf('../../Work/GIS/cb_2018_us_state_20m/or_wa_borders.shp') |> 
#   st_transform("EPSG:3741")
# template_general <- crop(template_general, or_wa)
# rm(or_wa)
makeTiles(template_general, ceiling(dim(template_general)/150)[1:2], 
          filename = paste0('data/tiles/',file_path_sans_ext(template_suffix),'_.tif'), 
          extend = TRUE, na.rm = TRUE, overwrite = TRUE)
rm(template_general)
# ------------------------------------------------------------------------------

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
# RUN THIS: 
plot(rast(file.path('data','tiles',"_topo_1981-2010_1301.tif")))
#list.files('data/tiles')
"_topo_1981-2010_1301.tif"|>
  walk(\(tile){
    
    tic('Per tile step')
    
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
 
    # Fine extract for all variables, at centroid, returns df with row for each
    # ----- This is failing for some reason... returning null/0. terra::extact works
    # since extract one point at a time with terra is fast, moving that below.
    # fine_dat <- exact_extract(template_fine, fine_centroids, include_xy = T, progress = T)

    # Preallocate dataframe
    # --------------- POPOULATE WITH X, Y ?
    result_df <- data.frame('pt_ID' = numeric(0), 'var' = character(0), 'gids' = numeric(0))
    
    toc()
    # --------------------------------------------------------------------------

    # Best approach appears to be a loop...! Through each element of coarse and template lists, 
    # which are dataframes, and each row of focal point centroids, which is a dataframe
    tic('Per point step')
    
    for (i in 1:length(coarse_dat)){ #141970
      
      # Distance calculations for each point in this tile
      x_dist <- st_coordinates(fine_centroids)[i,'X'] - coarse_dat[[i]][['x']]
      y_dist <- st_coordinates(fine_centroids)[i,'Y'] - coarse_dat[[i]][['y']]
      dists <- sqrt( x_dist^2 + y_dist^2 )
      coarse_dat[[i]][['di']] <- dists
      
      focal_result_df <- list('def','aet','tmin','tmax') |> # REPLACE HERE ALL THE YEARS... 
        map_dfr(\(clim_var){
          
          # Fine-grid information from fine-grid focal location
          Xp <- st_coordinates(fine_centroids[i,])[,'X']
          Yp <- st_coordinates(fine_centroids[i,])[,'Y']
          Cp <- terra::extract(template_fine, fine_centroids[i,])[[clim_var]]
          
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
          
          # LM version ---------------------------------------------------------
          # lm_mod <- lm(Zi ~ Xi + Yi + Ci, data = model_df) 
          
          x_lm <- model.matrix(Zi~Xi+Yi+Ci, data = model_df) 
          y_lm <- model_df[,'Zi']
          lm_mod <- .lm.fit(x_lm, y_lm)
          
          Cx <- lm_mod$coefficients[2]
          Cy <- lm_mod$coefficients[3]
          Cc <- lm_mod$coefficients[4]
          
          # Inverse distance weighting
          # This forumula is provided on page 5 of Flint and Flint 2012
          sum1 <- sum(( model_df$Zi + (Xp-model_df$Xi)*Cx + (Yp-model_df$Yi)*Cy + (Cp-model_df$Ci)*Cc ) / model_df$di^2)
          sum2 <- sum(1/model_df$di^2)
          Z = sum1/sum2 
          # --------------------------------------------------------------------
          
          # Return as dataframe
          return(data.frame('var' = clim_var,
                            'gids' = Z))
        })
      
      result_df <- rbind(result_df, focal_result_df)
    }
    
    # Pivot the result to wide, for easy write-out as raster stack
    result_df_wide <- result_df |> 
      pivot_wider(names_from = var, values_from = gids) 
    
    # Write out this tile as a tiff
    st_coordinates(fine_centroids) |> 
      cbind(result_df_wide) |> 
      rast(type = 'xyz', crs = tile_rast) |> 
      writeRaster(paste0('data/gids_output/','ds-270-m_',c('def','aet','tmin','tmax'),tile), 
                  overwrite = T)
    toc()
  })

plot(rast("data/gids_output/ds-270-m_def_topo_1981-2010_578.tif"))
    
    
            
    
              
              
              
           


toc()
plan(sequential)

# DELETE TILES FOR THIS clim_varIABLE!

# Testing info
complete <- rast(file.path(cogs_path,paste0('def',template_suffix)))
tile <- crop(complete, bounds)
(242.15   *(ncell(complete)/ncell(tile))) / 60 / 60

# Approximately 48 compute hours per clim_variable for entire western US  


