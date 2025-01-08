library(terra)
library(magrittr)
library(sf)
library(purrr)
# Time periods
times <- paste0("hist_", 1961:2022) # ,,'1961-1990',paste0('2C_',1985:2015)
# times <- paste0('2C_',1985:2015)
out_dir <- "merged_output"

if (!file.exists(paste0("data/", out_dir))) {
  dir.create(paste0("data/", out_dir))
}

crop_layer <- st_read("data/crop_layer_wgs.gpkg") |>
  st_transform(crs = "EPSG:4326")

times %>%
  walk(\(time){
    print(paste0("Starting ", time, "..."))

    file_rast <- grep(
      time,
      list.files(
        "data/gids_output",
        full.names = T
      ),
      value = T
    ) %>%
      map(., rast) %>%
      sprc() %>%
      merge()

    # AET
    file_rast[["aet"]][file_rast[["aet"]] < 0] <- 0
    file_rast[["aet"]][file_rast[["aet"]] > 1000] <- 1000
    file_rast[["aet"]] <- round(file_rast[["aet"]], 1)
    # DEF
    file_rast[["def"]][file_rast[["def"]] < 0] <- 0
    file_rast[["def"]][file_rast[["def"]] > 2500] <- 2500
    file_rast[["def"]] <- round(file_rast[["def"]], 1)
    # TMAX
    file_rast[["tmax"]][file_rast[["tmax"]] < -5] <- -5 # This happens to not be needed, but just for completeness
    file_rast[["tmax"]][file_rast[["tmax"]] > 50] <- 50
    file_rast[["tmax"]] <- round(file_rast[["tmax"]], 2)
    # TMIN
    file_rast[["tmin"]][file_rast[["tmin"]] < -25] <- -25
    file_rast[["tmin"]][file_rast[["tmin"]] > 20] <- 20
    file_rast[["tmin"]] <- round(file_rast[["tmin"]], 2)



    file_rast <- crop(x = file_rast, y = crop_layer, mask = TRUE)

    print(paste0("Writing ", time, "..."))

    writeRaster(file_rast,
      paste0("data/", out_dir, "/topoterra_", time, ".tif"),
      datatype = "FLT4S",
      gdal = c(
        "TILED=YES",
        "BLOCKXSIZE=128",
        "BLOCKYSIZE=128",
        "COMPRESS=DEFLATE"
      ),
      overwrite = TRUE
    )
    rm(file_rast)
    gc()
  })
