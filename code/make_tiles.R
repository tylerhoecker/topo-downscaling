library(terra)
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
