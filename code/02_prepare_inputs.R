# ------------------------------------------------------------------------------
# Prepared 2/27/2023 by Tyler Hoecker: https://github.com/tylerhoecker
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# This script prepares the geospatial layers needed to downscale future climate projections
# It performs the following steps:
# - Reproject raster grids from whatever they are in to ESPG 3741, UTM Zone 11 
# - Mask out areas that are known to be covered by water otherwise could never support vegetation
# - Reformat tiff as a Cloud-Optimized-Geotiff (COG), a preferred format for cloud-based computing
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load packages
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs,require,character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}
using('terra','gdalUtilities','purrr')

# ------------------------------------------------------------------------------
# Assuming you opened the script in the 'code' directory, set to root, or whatever relative path you prefer
# ------------------------------------------------------------------------------
# Climate variables 
files <- list.files(path = 'data/climate_inputs/', pattern = '.*tif$')

# Load the 'water mask' based on landfire, which was created and transformed to the same CRS
#lf_water_mask <- rast('data/parks_water_mask_250.tif')
lf_water_mask <- rast('data/lf_water_mask_4326.tiff')
# Resample water mask to an example of topo and terra so it matches each
grid_topo <- rast('data/climate_inputs/tmin_topo_hist_1981-2010.tif')
# grid_topo_utm <- project(grid_topo, "EPSG:32611")
# lf_water_mask_topo <- resample(lf_water_mask, grid_topo_utm)

# # Resample water mask to an example of topo and terra so it matches each
# grid_terra <- rast('../data/climate_inputs/tmin_terra_1961-1990.tif')
# grid_terra_utm <- project(grid_terra, "EPSG:32611")
# lf_water_mask_terra <- resample(lf_water_mask, grid_terra_utm)
lf_water_mask <- resample(lf_water_mask, grid_topo, "near")

# Do the step outlined above for every file in the climate_inputs directory
files |> 
  walk(function(file){
    
    # Original file format
    grid_original <- rast(paste0('data/climate_inputs/',file))
    
    # Make sure all share the same projection
    # grid_utm <- project(grid_original, "EPSG:32611")
    
    # Mask out water, if TopoFire data
    dataset <- sub('_[0-9]+.*','', sub('[a-z]+_','',file))
    
    if(dataset == 'topo_hist'){
      grid_final <- mask(grid_original, lf_water_mask)
    }else{
      grid_final <- crop(grid_original, grid_topo)
    }

    # Write it out
    writeRaster(grid_final, paste0('data/cogs/',file),
                overwrite = TRUE, 
                gdal = c("TILED=YES",
                         "COPY_SRC_OVERVIEWS=YES",
                         "COMPRESS=DEFLATE"))
   
   
    
  })


# #-------------------------------------------------------------------------------
# # Make the fine point template COG (all points w/o water mask that will be used to create tiles)
# #-------------------------------------------------------------------------------
# # Original file format
# file <- 'def_topo_1981-2010.tif'
# 
# grid_original <- rast(file.path('../data/cogs/',file))
# 
# # Make sure all share the same projection
# #grid_original <- project(grid_original, grid_master)
# grid_original_utm <- project(grid_original, "EPSG:3741")
# 
# # Turn all values to 1 to avoid confusion with real data
# grid_original_utm[!is.na(grid_original_utm)] <- 1
# 
# # Not cropping because data from Canada and Mexico can be used
# # Write it out
# writeRaster(grid_original_utm, paste0('temp_',file),
#             overwrite = TRUE)
# 
# gdal_translate(
#   src_dataset = paste0('temp_',file),
#   # This is different
#   dst_dataset = '../data/master_template.tif',
#   co = matrix(c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=DEFLATE"),ncol = 1)
# )
# 
# unlink(paste0('temp_',file))






