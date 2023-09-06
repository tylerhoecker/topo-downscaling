library(terra)
library(tools)
library(terra)
library(purrr)
# Layer names
# These are both historical datasets. 1st step is to downscale historical TerraClimate
template_suffix <- '_topo_1981-2010.tif' 
coarse_suffix <- '_terra_1981-2010.tif'

# Make tiles -------------------------------------------------------------------
# Build information about fine-scale coordinates for this tile
# Not necessary to build tile for coarse raster
# Climate variable for this step is arbitrary - spatial info identical for all
template_general <- rast(file.path('data','cogs',paste0('def',template_suffix)))

makeTiles(template_general, ceiling(dim(template_general)/100)[1:2],
          filename = paste0('data/tiles/',file_path_sans_ext(template_suffix),'_.tif'),
          extend = TRUE, na.rm = TRUE, overwrite = TRUE)
rm(template_general)

# Eventually, do this for each year...
clim_vars <- c('def','aet','tmin','tmax')

# Inputs for entire western US, do once
# Load a complete coarse template 
coarse_files <- paste0(clim_vars, coarse_suffix) #,'aet','tmin','tmax'
ds_coarse <- rast(file.path('data','cogs', coarse_files))
names(ds_coarse) <- clim_vars 

# Create a coarse (4 km) version of the "template", which will align with future data,
# to regress against the coarse data
fine_files <- paste0(clim_vars, template_suffix) #,'aet','tmin','tmax'
template_fine <- rast(file.path('data','cogs', fine_files))
names(template_fine) <- clim_vars #,

# Crop Terra (global) to Topo (N. America) - will still include parts of Mexico and CA
ds_coarse <- crop(ds_coarse, template_fine)

# Resample, as a form of aggregation, the fine historical data to the desired coarse grid
template_coarse <- resample(template_fine, ds_coarse, 'bilinear')
# Not sure why this is necessary... but it is
crs(template_coarse) <- crs(ds_coarse)
names(template_coarse) <- paste0(clim_vars)

# Crop data layers to area around each tile and save...
list.files('data/tiles') |> 
  walk(\(tile){
    print(paste0('Running:',tile))
    tile_rast <-  rast(paste0('data/tiles/',tile))
    tile_border <- ext(tile_rast)
    tile_border <- vect(tile_border)  
    tile_border <- buffer(tile_border, width = 50000)
    ds_coarse_tile <- crop(ds_coarse, tile_border)
    template_fine_tile <- crop(template_fine, tile_border)
    template_coarse_tile <- crop(template_coarse, tile_border)
    writeRaster(ds_coarse_tile, paste0('data/tile_templates/ds_coarse',tile), overwrite = T)
    writeRaster(template_fine_tile, paste0('data/tile_templates/template_fine',tile), overwrite = T)
    writeRaster(template_coarse_tile, paste0('data/tile_templates/template_coarse',tile), overwrite = T)
  })




