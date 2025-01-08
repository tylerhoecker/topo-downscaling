library(terra)
library(tidyverse)

### Historical Normal
reference_crs <- crs(rast("data/merged_output/topoterra_hist_1961.tif"))
#### TerraClimate
years <- 1961:2022
terraclimate_files_hist <- list.files("data/cogs", pattern = paste0("terra_hist_[0-9]{4}.tif$"), full.names = TRUE)
aet <- terraclimate_files_hist[grep("aet", terraclimate_files_hist)] %>%
  rast() %>%
  mean()
def <- terraclimate_files_hist[grep("def", terraclimate_files_hist)] %>%
  rast() %>%
  mean()
tmax <- terraclimate_files_hist[grep("tmax", terraclimate_files_hist)] %>%
  rast() %>%
  mean()
tmin <- terraclimate_files_hist[grep("tmin", terraclimate_files_hist)] %>%
  rast() %>%
  mean()

terraclimate_normal <- c(aet, def, tmax, tmin) %>%
  project(reference_crs)
names(terraclimate_normal) <- c("aet", "def", "tmax", "tmin")

##  clean terraclim like topoterra
# AET
terraclimate_normal[["aet"]][terraclimate_normal[["aet"]] < 0] <- 0
terraclimate_normal[["aet"]][terraclimate_normal[["aet"]] > 1000] <- 1000
terraclimate_normal[["aet"]] <- round(terraclimate_normal[["aet"]], 1)
# DEF
terraclimate_normal[["def"]][terraclimate_normal[["def"]] < 0] <- 0
terraclimate_normal[["def"]][terraclimate_normal[["def"]] > 2500] <- 2500
terraclimate_normal[["def"]] <- round(terraclimate_normal[["def"]], 1)
# TMAX
terraclimate_normal[["tmax"]][terraclimate_normal[["tmax"]] < -5] <- -5 # This happens to not be needed, but just for completeness
terraclimate_normal[["tmax"]][terraclimate_normal[["tmax"]] > 50] <- 50
terraclimate_normal[["tmax"]] <- round(terraclimate_normal[["tmax"]], 2)
# TMIN
terraclimate_normal[["tmin"]][terraclimate_normal[["tmin"]] < -25] <- -25
terraclimate_normal[["tmin"]][terraclimate_normal[["tmin"]] > 20] <- 20
terraclimate_normal[["tmin"]] <- round(terraclimate_normal[["tmin"]], 2)

# write out
writeRaster(terraclimate_normal, "data/merged_output/terra_hist_1961-2022.tif",
  overwrite = TRUE,
  datatype = "FLT4S"
)

#### TopoTerra
# load topoterra
topoterra_files_hist <- list.files("data/merged_output", pattern = "topoterra_hist_[0-9]{4}.tif", full.names = TRUE)
# filter to 1961:1:1990
topoterra_files_hist <- topoterra_files_hist
topoterra_stack <- map(topoterra_files_hist, rast)
# organize into indivudal stacks for each variable
aet <- vector("list", length(topoterra_stack))
def <- vector("list", length(topoterra_stack))
tmax <- vector("list", length(topoterra_stack))
tmin <- vector("list", length(topoterra_stack))
for (i in 1:length(topoterra_stack)) {
  aet[[i]] <- topoterra_stack[[i]][[1]]
  def[[i]] <- topoterra_stack[[i]][[2]]
  tmax[[i]] <- topoterra_stack[[i]][[3]]
  tmin[[i]] <- topoterra_stack[[i]][[4]]
}
# calculate normal all
aet_normal <- rast(aet) %>%
  mean() %>%
  round(1)
def_normal <- rast(def) %>%
  mean() %>%
  round(1)
tmax_normal <- rast(tmax) %>%
  mean() %>%
  round(2)
tmin_normal <- rast(tmin) %>%
  mean() %>%
  round(2)
# calculate normal 1961:1990
aet_normal_1961_1990 <- rast(aet[1:30]) %>%
  mean() %>%
  round(1)
def_normal_1961_1990 <- rast(def[1:30]) %>%
  mean() %>%
  round(1)
tmax_normal_1961_1990 <- rast(tmax[1:30]) %>%
  mean() %>%
  round(2)
tmin_normal_1961_1990 <- rast(tmin[1:30]) %>%
  mean() %>%
  round(2)



# combine
topoterra_normal <- c(aet_normal, def_normal, tmax_normal, tmin_normal)
names(topoterra_normal) <- names(topoterra_stack[[1]])
topoterra_normal_1961_1990 <- c(aet_normal_1961_1990, def_normal_1961_1990, tmax_normal_1961_1990, tmin_normal_1961_1990)
names(topoterra_normal_1961_1990) <- names(topoterra_stack[[1]])
# write out
writeRaster(topoterra_normal, "data/merged_output/topoterra_hist_1961-2022.tif",
  overwrite = TRUE,
  gdal = c(
    "TILED=YES",
    "COPY_SRC_OVERVIEWS=YES",
    "COMPRESS=DEFLATE"
  )
)
writeRaster(topoterra_normal_1961_1990, "data/merged_output/topoterra_hist_1961-1990.tif",
  overwrite = TRUE,
  gdal = c(
    "TILED=YES",
    "COPY_SRC_OVERVIEWS=YES",
    "COMPRESS=DEFLATE"
  )
)

# make 1985-2015 normal of hist topoterra
years <- 1985:2015
topoterra_files_hist <- list.files("data/merged_output", pattern = "topoterra_hist_[0-9]{4}.tif", full.names = TRUE)
topoterra_files_hist <- topoterra_files_hist[grepl(paste0(years, collapse = "|"), topoterra_files_hist)]
topoterra_stack <- map(topoterra_files_hist, rast)
# organize into indivudal stacks for each variable
aet <- vector("list", length(topoterra_stack))
def <- vector("list", length(topoterra_stack))
tmax <- vector("list", length(topoterra_stack))
tmin <- vector("list", length(topoterra_stack))
for (i in 1:length(topoterra_stack)) {
  aet[[i]] <- topoterra_stack[[i]][[1]]
  def[[i]] <- topoterra_stack[[i]][[2]]
  tmax[[i]] <- topoterra_stack[[i]][[3]]
  tmin[[i]] <- topoterra_stack[[i]][[4]]
}
# calculate normal
aet_normal <- rast(aet) %>%
  mean() %>%
  round(1)
def_normal <- rast(def) %>%
  mean() %>%
  round(1)
tmax_normal <- rast(tmax) %>%
  mean() %>%
  round(2)
tmin_normal <- rast(tmin) %>%
  mean() %>%
  round(2)

# combine
topoterra_normal <- c(aet_normal, def_normal, tmax_normal, tmin_normal)
names(topoterra_normal) <- c("aet", "def", "tmax", "tmin")

# write out
writeRaster(topoterra_normal, "data/merged_output/topoterra_hist_1985-2015.tif",
  overwrite = TRUE,
  gdal = c(
    "TILED=YES",
    "COPY_SRC_OVERVIEWS=YES",
    "COMPRESS=DEFLATE"
  )
)


### 2C Normal

#### TerraClimate
terraclimate_files_2c <- list.files("data/cogs", pattern = "2C_[0-9]{4}.tif$", full.names = T)

aet <- terraclimate_files_2c[grep("aet", terraclimate_files_2c)] %>%
  rast() %>%
  mean()
def <- terraclimate_files_2c[grep("def", terraclimate_files_2c)] %>%
  rast() %>%
  mean()
tmax <- terraclimate_files_2c[grep("tmax", terraclimate_files_2c)] %>%
  rast() %>%
  mean()
tmin <- terraclimate_files_2c[grep("tmin", terraclimate_files_2c)] %>%
  rast() %>%
  mean()

terraclimate_normal <- c(aet, def, tmax, tmin) %>%
  project(reference_crs)
names(terraclimate_normal) <- c("aet", "def", "tmax", "tmin")

##  clean terraclim like topoterra
# AET
terraclimate_normal[["aet"]][terraclimate_normal[["aet"]] < 0] <- 0
terraclimate_normal[["aet"]][terraclimate_normal[["aet"]] > 1000] <- 1000
terraclimate_normal[["aet"]] <- round(terraclimate_normal[["aet"]], 1)
# DEF
terraclimate_normal[["def"]][terraclimate_normal[["def"]] < 0] <- 0
terraclimate_normal[["def"]][terraclimate_normal[["def"]] > 2500] <- 2500
terraclimate_normal[["def"]] <- round(terraclimate_normal[["def"]], 1)
# TMAX
terraclimate_normal[["tmax"]][terraclimate_normal[["tmax"]] < -5] <- -5 # This happens to not be needed, but just for completeness
terraclimate_normal[["tmax"]][terraclimate_normal[["tmax"]] > 50] <- 50
terraclimate_normal[["tmax"]] <- round(terraclimate_normal[["tmax"]], 2)
# TMIN
terraclimate_normal[["tmin"]][terraclimate_normal[["tmin"]] < -25] <- -25
terraclimate_normal[["tmin"]][terraclimate_normal[["tmin"]] > 20] <- 20
terraclimate_normal[["tmin"]] <- round(terraclimate_normal[["tmin"]], 2)

# write out
writeRaster(terraclimate_normal, "data/merged_output/terra_2C_1985-2015.tif",
  overwrite = TRUE,
  datatype = "FLT4S"
)

#### TopoTerra
# load topoterra
topoterra_files_2c <- list.files("data/merged_output", pattern = "topoterra_2C_[0-9]{4}.tif", full.names = TRUE)
topoterra_stack <- map(topoterra_files_2c, rast)
# organize into indivudal stacks for each variable
aet <- vector("list", length(topoterra_stack))
def <- vector("list", length(topoterra_stack))
tmax <- vector("list", length(topoterra_stack))
tmin <- vector("list", length(topoterra_stack))
for (i in 1:length(topoterra_stack)) {
  aet[[i]] <- topoterra_stack[[i]][[1]]
  def[[i]] <- topoterra_stack[[i]][[2]]
  tmax[[i]] <- topoterra_stack[[i]][[3]]
  tmin[[i]] <- topoterra_stack[[i]][[4]]
}
# calculate normal
aet_normal <- rast(aet) %>%
  mean() %>%
  round(1)
def_normal <- rast(def) %>%
  mean() %>%
  round(1)
tmax_normal <- rast(tmax) %>%
  mean() %>%
  round(2)
tmin_normal <- rast(tmin) %>%
  mean() %>%
  round(2)

# combine
topoterra_normal <- c(aet_normal, def_normal, tmax_normal, tmin_normal)
names(topoterra_normal) <- names(topoterra_stack[[1]])

# write out
writeRaster(topoterra_normal, "data/merged_output/topoterra_2C_1985-2015.tif",
  overwrite = TRUE,
  gdal = c(
    "TILED=YES",
    "COPY_SRC_OVERVIEWS=YES",
    "COMPRESS=DEFLATE"
  )
)

# make delta normal
topoterra_2c <- rast("data/merged_output/TopoTerra_2C_1985-2015.tif")
topoterra_hist <- rast("data/merged_output/TopoTerra_1961-2022.tif")
topoterra_delta <- topoterra_2c - topoterra_hist
# round to the same as above
topoterra_delta[["aet"]] <- round(topoterra_delta[["aet"]], 1)
topoterra_delta[["def"]] <- round(topoterra_delta[["def"]], 1)
topoterra_delta[["tmax"]] <- round(topoterra_delta[["tmax"]], 2)
topoterra_delta[["tmin"]] <- round(topoterra_delta[["tmin"]], 2)

writeRaster(topoterra_delta, "data/merged_output/TopoTerra_delta_all_years.tif",
  overwrite = TRUE,
  gdal = c(
    "TILED=YES",
    "COPY_SRC_OVERVIEWS=YES",
    "COMPRESS=DEFLATE"
  )
)

# make delta normal of topoterra over same tinme period 1985:2015
topoterra_hist <- rast("data/merged_output/TopoTerra_1985-2015.tif")
topoterra_2c <- rast("data/merged_output/TopoTerra_2C_1985-2015.tif")
topoterra_delta <- topoterra_2c - topoterra_hist
# round to the same as above
topoterra_delta[["aet"]] <- round(topoterra_delta[["aet"]], 1)
topoterra_delta[["def"]] <- round(topoterra_delta[["def"]], 1)
topoterra_delta[["tmax"]] <- round(topoterra_delta[["tmax"]], 2)
topoterra_delta[["tmin"]] <- round(topoterra_delta[["tmin"]], 2)

writeRaster(topoterra_delta, "data/merged_output/TopoTerra_delta_1985-2015.tif",
  overwrite = TRUE,
  gdal = c(
    "TILED=YES",
    "COPY_SRC_OVERVIEWS=YES",
    "COMPRESS=DEFLATE"
  )
)

# delta normal of topoterra from hist_1961-1990 to 2C 1985-2015
topoterra_hist <- rast("data/merged_output/topoterra_hist_1961-1990.tif")
topoterra_2c <- rast("data/merged_output/topoterra_2C_1985-2015.tif")
topoterra_delta <- topoterra_2c - topoterra_hist
# round to the same as above
topoterra_delta[["aet"]] <- round(topoterra_delta[["aet"]], 1)
topoterra_delta[["def"]] <- round(topoterra_delta[["def"]], 1)
topoterra_delta[["tmax"]] <- round(topoterra_delta[["tmax"]], 2)
topoterra_delta[["tmin"]] <- round(topoterra_delta[["tmin"]], 2)

writeRaster(topoterra_delta, "data/merged_output/TopoTerra_delta_1961-1990_1985-2015.tif",
  overwrite = TRUE,
  gdal = c(
    "TILED=YES",
    "COPY_SRC_OVERVIEWS=YES",
    "COMPRESS=DEFLATE"
  )
)



## quick test. go through every file in merged_output and check if
# aet and def are rounded to 1 decimal place and tmax and tmin are rounded to 2 decimal places
# if not, flag them
merged_output_files <- list.files("data/merged_output", pattern = "topoterra_.*tif$", full.names = TRUE)
for (file in merged_output_files) {
  r <- rast(file)

  suppressWarnings({
    if (isTRUE(all.equal(round(r[["aet"]], 1), r[["aet"]]))) {
      print(paste0(file, " aet"))
    }
    if (isTRUE(all.equal(round(r[["def"]], 1), r[["def"]]))) {
      print(paste0(file, " def"))
    }
    if (isTRUE(all.equal(round(r[["tmax"]], 2), r[["tmax"]]))) {
      print(paste0(file, " tmax"))
    }
    if (isTRUE(all.equal(round(r[["tmin"]], 2), r[["tmin"]]))) {
      print(paste0(file, " tmin"))
    }
  })
  rm(r)
  gc()
}


# rename based on the new naming convention
## rules topoterra -> TopoTerra
## hist is removed from the name, but 2C is kept
## TopoTerra_2C_year.tif
## keep the range for normals
## ex. TopoTerra_1961-2022.tif

all_topoterra_files <- list.files("data/merged_output", pattern = "topoterra_.*tif$", full.names = TRUE)
all_topoterra_files %>%
  walk(\(x) {
    new_name <- x %>%
      str_replace_all("topoterra", "TopoTerra") %>%
      str_remove("hist_")
    file.rename(x, new_name)
  })

terra_names <- list.files("data/merged_output", pattern = "terra_.*tif$", full.names = TRUE)
terra_names %>%
  walk(\(x) {
    new_name <- x %>%
      str_replace_all("terra", "Terra") %>%
      str_remove("hist_")
    file.rename(x, new_name)
  })
