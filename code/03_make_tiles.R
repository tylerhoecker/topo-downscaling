# This script breaks the input datasets into tiles, so that the downscaling procedure
# can be run in parallel on multiple cores.
# This script also does some cropping so that data outside the study area are not
# read during the downscaling procedure, which uses RAM.

# Load packages
library(terra)
library(purrr)

clim_vars <- c('aet','def','tmax','tmin')

# A 'suffix' or naming convention common to all files to be downscaled
coarse_name <- 'terra'

# Time periods
times <- c('1961-1990','2C_1985-2015') #,,paste0('2C_',1985:2015)

# A 'suffix' or naming convention common to the files be used as a 'template' for downscaling
templ_name <- '_topo_1981-2010.tif'

# Buffer distance for geospatial regression, in meters
buff_dist <- 50000


# Make tiles -------------------------------------------------------------------
# Build information about fine-scale coordinates for this tile
# Climate variable for this step is arbitrary - spatial info identical for all
template_general <- rast(paste0('../data/cogs/',clim_vars[1],templ_name))

dir.create('../data/tiles')

makeTiles(template_general, ceiling(dim(template_general)/200)[1:2],
          filename = paste0('../data/tiles/_.tif'),
          extend = TRUE, na.rm = TRUE, overwrite = TRUE)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Crop data layers to area around each tile and save...
#-------------------------------------------------------------------------------

# Create a output directory for these
if (!file.exists('../data/tile_templates')) {
  dir.create('../data/tile_templates')
}

as.list(times) %>% 
  walk(\(time){
    
    print(paste0('Running:',time))
    
    # Stack of coarse files (all variables) for this time --
    ds_coarse <- rast(paste0('../data/cogs/',clim_vars,'_',coarse_name,'_',time,'.tif'))
    names(ds_coarse) <- clim_vars
    ds_coarse <- round(ds_coarse, 1)
    
    # Fine template 
    template_fine <- rast(paste0('../data/cogs/',clim_vars,templ_name))
    names(template_fine) <- clim_vars
    template_fine <- round(template_fine, 1)
    
    # Resample, as a form of aggregation, the fine data to coarse grid, so they can be regressed
    template_coarse <- resample(template_fine, ds_coarse, 'bilinear', threads = T)
    
    
    list.files('../data/tiles') %>% 
      walk(\(tile){ # tile=list.files('../data/tiles')[[100]]
        
        # Progress message
        print(paste0('Running:',tile))
        
        # Load the tile
        tile_rast <- rast(paste0('../data/tiles/',tile))
        
        # Turn the extent of the raster into a spatial vector
        tile_border <- ext(tile_rast)
        tile_border <- vect(tile_border)
        
        # Create a buffer around the tile (in this case, 50 km )
        tile_border <- buffer(tile_border, width = buff_dist)
        
        # Crop input rasters (all variable layers at once)
        ds_coarse_tile <- crop(ds_coarse, tile_border)
        template_fine_tile <- crop(template_fine, tile_border)
        template_coarse_tile <- crop(template_coarse, tile_border)
        
        # Write them out
        writeRaster(ds_coarse_tile, 
                    paste0('../data/tile_templates/ds_coarse_',time,tile))
        writeRaster(template_fine_tile, 
                    paste0('../data/tile_templates/template_fine_',time,tile))
        writeRaster(template_coarse_tile, 
                    paste0('../data/tile_templates/template_coarse_',time,tile))
        
        rm(tile_rast, tile_border, ds_coarse_tile, template_fine_tile, template_coarse_tile)
        gc()
      })
    rm(ds_coarse, template_fine, template_coarse)
    gc()
  })





