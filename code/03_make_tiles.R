# This script breaks the input datasets into tiles, so that the downscaling procedure
# can be run in parallel on multiple cores.
# This script also does some cropping so that data outside the study area are not
# read during the downscaling procedure, which uses RAM.

# Load packages
library(terra)
#library(tools)
library(terra)
library(purrr)

# Make tiles -------------------------------------------------------------------
# Build information about fine-scale coordinates for this tile
# Climate variable for this step is arbitrary - spatial info identical for all

# NEW - Use a version of topoclimate to build the template, so that we predict over 
# all points, but don't use non-veg points in regressions
template_general <- rast('../data/cogs/def_topo_1981-2010.tif')

dir.create('../data/tiles')

makeTiles(template_general, ceiling(dim(template_general)/50)[1:2],
          filename = paste0('../data/tiles/_.tif'),
          extend = TRUE, na.rm = TRUE, overwrite = TRUE)
rm(template_general)
#-------------------------------------------------------------------------------

# Crop data layers to area around each tile and save...
# dir.create('data/tile_templates')
# 
# 
# list.files('../data/tiles') |> 
#   walk(\(tile){
#     print(paste0('Running:',tile))
#     tile_rast <-  rast(paste0('data/tiles/',tile))
#     tile_border <- ext(tile_rast)
#     tile_border <- vect(tile_border)  
#     tile_border <- buffer(tile_border, width = 50000)
#     ds_coarse_tile <- crop(ds_coarse, tile_border)
#     template_fine_tile <- crop(template_fine, tile_border)
#     template_coarse_tile <- crop(template_coarse, tile_border)
#     writeRaster(ds_coarse_tile, paste0('data/tile_templates/ds_coarse',tile), overwrite = T)
#     writeRaster(template_fine_tile, paste0('data/tile_templates/template_fine',tile), overwrite = T)
#     writeRaster(template_coarse_tile, paste0('data/tile_templates/template_coarse',tile), overwrite = T)
#   })




