library(terra)
library(tidyverse)

# Time periods
times <- paste0(1968:1990) #,,'1961-1990',paste0('2C_',1985:2015)

out_dir <- 'merged_output'

if (!file.exists(paste0('data/',out_dir))) {
  dir.create(paste0('data/',out_dir))
}

crop_layer <- st_read("data/crop_layer_wgs.gpkg") |>
    st_transform(crs = "EPSG:3857")

as.list(times) %>%
  walk(\(time){

    print(paste0("Starting ", time, "..."))

    file_rast <- grep(
      time,
      list.files(
        'data/gids_output', 
        full.names = T
      ),
      value = T
      ) %>% 
      map(., rast) %>% 
      sprc() %>% 
      merge()

    # AET
    file_rast[['aet']][file_rast[['aet']] < 0] = 0
    file_rast[['aet']][file_rast[['aet']] > 1000] = 1000
    file_rast[['aet']] <- round(file_rast[['aet']],0) 
    # DEF
    file_rast[['def']][file_rast[['def']] < 0] = 0
    file_rast[['def']][file_rast[['def']] > 2500] = 2500
    file_rast[['def']] <- round(file_rast[['def']],0)
    # TMAX
    file_rast[['tmax']][file_rast[['tmax']] < -5] = -5 # This happens to not be needed, but just for completeness
    file_rast[['tmax']][file_rast[['tmax']] > 35] = 35
    file_rast[['tmax']] <- round(file_rast[['tmax']], 2)*10^2
    # TMIN
    file_rast[['tmin']][file_rast[['tmin']] < -25] = -25
    file_rast[['tmin']][file_rast[['tmin']] > 20] = 20
    file_rast[['tmin']] <- round(file_rast[['tmin']], 2)*10^2

  file_rast <- file_rast |>
    terra::project("EPSG:3857") #4326 works

  file_rast <- crop(x = file_rast, y = crop_layer, mask = TRUE)

  print(paste0("Writing ", time, "..."))

  writeRaster(file_rast, 
              paste0("data/", out_dir, "/topoterra_hist_", time, '.tif'), 
              datatype = "INT2S",
              gdal = c("PROJECTION=EPSG:3857",
                       "TILED=YES",
                       "BLOCKXSIZE=128",
                       "BLOCKYSIZE=128",
                       "OVERVIEW-RESAMPLING=NEAREST",
                       "COMPRESS=DEFLATE"),
              overwrite = TRUE)
    

  })


