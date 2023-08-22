# ------------------------------------------------------------------------------
# Prepared 2/27/2023 by Tyler Hoecker: https://github.com/tylerhoecker
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# This script prepares the geospatial layers needed to downscale future climate projections
#
# 250-m climate data from TopoFire (Holden et al.) will be used to downscale
# 4-km projected data, which was downscaled and bias-corrected from GCM output 
# by TerraClimate (Abatzoglou et al.)
# The downscaling process will be conducted in a cloud computing environment (CyVerse).

# The required files, all in cloud-optimized geotiff (COG) format, are:
# - Fine-scale "template" (in this case, 250-m climate data from TopoFire [Holden et al.])
# - Coarse-scale "template" that will be regressed against coarse-scale target layer
# - Coarse-scale "target" layer that will be downscaled
# - Multiply above by the number of climate variables / target layers to be downscaled.
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

# ------------------------------------------------------------------------------
# Data paths (may not generalize)
# Path where the tifs in non-cog format are located
# Do this for 'terraclimate' and 'topofire' rasters
file_path <- file.path('data','terraclimate')
# Path where you want to save tifs in COG format
cogs_path <- file.path('data','cogs')

# Climate variables 
files <- list.files(path = file_path, pattern = '.*tif$')


# MAYBE NOT...
# # Master grid - decide what all grids should be cropped to, so all extents match
# grid_master <- rast('data/topofire/aet_topo_1981-2010.tif')

files |> 
  map(function(file){
    # Original file format
    grid_original <- rast(file.path(file_path,file))
    
    # Make sure all share the same projection
    #grid_original <- project(grid_original, grid_master)
    grid_original_utm <- project(grid_original, "EPSG:3741")
   
    ## Crop the the master grid
    #grid_original_crop <- crop(grid_original, grid_master)
    writeRaster(grid_original_utm, paste0('temp_',file),
                overwrite = TRUE)
   
    gdal_translate(
      src_dataset = paste0('temp_',file),
      dst_dataset = file.path(cogs_path,file),
      co = matrix(c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=DEFLATE"),ncol = 1)
    )
   
    unlink(paste0('temp_',file))
  })



  





