# This script breaks the input datasets into tiles, so that the downscaling procedure
# can be run in parallel on multiple cores.
# This script also does some cropping so that data outside the study area are not
# read during the downscaling procedure, which uses RAM.

# Load packages
library(terra)
library(purrr)
library(furrr)

clim_vars <- c("aet","def","tmax","tmin")#c('aet','def','tmax','tmin')

# A 'suffix' or naming convention common to all files to be downscaled
coarse_name <- 'terra_2C'
coarse_name <- 'terra_hist'
# Time periods
times <- c(1985:2015) #paste0(1961:1990) #c('1961-1990','2C_1985-2015') #,,paste0('2C_',1985:2015)'1961-1990','2C_1985-2015'
times <- c(1961:2022)
# A 'suffix' or naming convention common to the files be used as a 'template' for downscaling
templ_name <- '_topo_hist_1981-2010.tif'

# Buffer distance for geospatial regression, in meters
buff_dist <- 50000

# Directory for tiles, relative to locatoin of this script
tile_dir <- 'data/tiles/'

# Directory for tile_templates (cropped rasters for each tile)
tile_templ_dir <- 'data/tile_templates/'

# Make tiles -------------------------------------------------------------------
#Build information about fine-scale coordinates for this tile
# Climate variable for this step is arbitrary - spatial info identical for all
template_general <- rast(paste0('data/cogs/',clim_vars[1],templ_name))

# Reduce area for testing - crop to PNW Ecoregions
pnw_ecoregions <- sf::st_read("data/pnw_ecoregions/pnw_ecoregions.gpkg")
template_general <- crop(template_general, pnw_ecoregions)
template_general <- mask(template_general, pnw_ecoregions)

dir.create(tile_dir)

makeTiles(template_general, ceiling(dim(template_general)/100)[1:2],
         filename = paste0(tile_dir,'_.tif'),
         extend = TRUE, na.rm = TRUE, overwrite = TRUE)
#-------------------------------------------------------------------------------


# write Raster function
writeCOG <- function(x, filename){
  writeRaster(x, 
              filename, 
              datatype = "INT2S",
              gdal = c("PROJECTION=EPSG:4326",
                       "TILED=YES",
                       "BLOCKXSIZE=128",
                       "BLOCKYSIZE=128",
                       "OVERVIEW-RESAMPLING=NEAREST",
                       "COMPRESS=DEFLATE"),
              overwrite = TRUE)
}

#-------------------------------------------------------------------------------
# Crop data layers to area around each tile and save...
#-------------------------------------------------------------------------------
# Create a output directory for these
if (!file.exists(tile_templ_dir)) {
  dir.create(tile_templ_dir)
}
 # Fine template 
    template_fine <- rast(paste0('data/cogs/',clim_vars,templ_name))
    names(template_fine) <- clim_vars
    template_fine <- round(template_fine, 1)
    
    # Resample, as a form of aggregation, the fine data to coarse grid, so they can be regressed
ds_coarse_template <- rast(paste0('data/cogs/',clim_vars,'_',coarse_name,'_',times[1],'.tif')) %>%
  project(crs(template_fine))
names(ds_coarse_template) <- clim_vars
    template_coarse <- resample(template_fine, ds_coarse_template, 'bilinear', threads = T)

#write out templates
    writeRaster(template_fine, paste0(tile_templ_dir,'template_fine_','.tif'), overwrite = TRUE)
    writeRaster(template_coarse, paste0(tile_templ_dir,'template_coarse_','.tif'), overwrite = TRUE)
        
plan(multicore, workers = 4)

times %>% 
  future_walk(\(time){
    
    n_complete <- length(
      list.files(
        tile_templ_dir, pattern = paste0('ds_coarse_', time, "_.*tif")
      )
    )

    n_tiles <- length(list.files(tile_dir))

    if( n_complete == n_tiles ){
      return(NULL)
    } else {
      print(paste0('Running:',time))
    }
    
    # Stack of coarse files (all variables) for this time --
    ds_coarse <- rast(paste0('data/cogs/',clim_vars,'_',coarse_name,'_',time,'.tif')) %>%
      project(crs("EPSG:4326"))
    names(ds_coarse) <- clim_vars
    ds_coarse <- round(ds_coarse, 1)
    
    #write out coarse data
    writeCOG(ds_coarse, paste0(tile_templ_dir,'ds_coarse_',coarse_name,"_",time,'.tif'))
   
    
    rm(ds_coarse, template_fine, template_coarse)
    gc()
  })





