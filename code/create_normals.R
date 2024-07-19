library(terra)
library(tidyverse)

### Historical Normal
reference_crs <- crs(rast("data/merged_output/topoterra_hist_1961.tif"))
#### TerraClimate
years <- 1961:2022
terraclimate_files_hist <- list.files("data/cogs", pattern = paste0("terra_hist_[0-9]{4}.tif$"), full.names = TRUE)
aet <- terraclimate_files_hist[grep("aet",terraclimate_files_hist)] %>% 
  rast() %>%
  median()
def <- terraclimate_files_hist[grep("def",terraclimate_files_hist)] %>%
  rast() %>%
  median()
tmax <- terraclimate_files_hist[grep("tmax",terraclimate_files_hist)] %>%
  rast() %>%
  median()
tmin <- terraclimate_files_hist[grep("tmin",terraclimate_files_hist)] %>%
  rast() %>%
  median()

terraclimate_normal <- c(aet, def, tmax, tmin) %>% 
  project(reference_crs)
names(terraclimate_normal) <- c("aet", "def", "tmax", "tmin")

##  clean terraclim like topoterra
# AET
terraclimate_normal[['aet']][terraclimate_normal[['aet']] < 0] = 0
terraclimate_normal[['aet']][terraclimate_normal[['aet']] > 1000] = 1000
terraclimate_normal[['aet']] <- round(terraclimate_normal[['aet']],0) 
# DEF
terraclimate_normal[['def']][terraclimate_normal[['def']] < 0] = 0
terraclimate_normal[['def']][terraclimate_normal[['def']] > 2500] = 2500
terraclimate_normal[['def']] <- round(terraclimate_normal[['def']],0)
# TMAX
terraclimate_normal[['tmax']][terraclimate_normal[['tmax']] < -5] = -5 # This happens to not be needed, but just for completeness
terraclimate_normal[['tmax']][terraclimate_normal[['tmax']] > 35] = 35
terraclimate_normal[['tmax']] <- round(terraclimate_normal[['tmax']], 2)*10^2
# TMIN
terraclimate_normal[['tmin']][terraclimate_normal[['tmin']] < -25] = -25
terraclimate_normal[['tmin']][terraclimate_normal[['tmin']] > 20] = 20
terraclimate_normal[['tmin']] <- round(terraclimate_normal[['tmin']], 2)*10^2

#write out
writeRaster(terraclimate_normal, "data/merged_output/terra_hist_1961-2022.tif",
            overwrite = TRUE, 
            gdal = c("TILED=YES",
                     "COMPRESS=DEFLATE"))

#### TopoTerra
#load topoterra
topoterra_files_hist <- list.files("data/merged_output", pattern = "topoterra_hist", full.names = TRUE) 
#filter to 1961:1:1990
topoterra_files_hist <- topoterra_files_hist
topoterra_stack <- map(topoterra_files_hist, rast)
#organize into indivudal stacks for each variable
aet <- vector("list", length(topoterra_stack))
def <- vector("list", length(topoterra_stack))
tmax <- vector("list", length(topoterra_stack))
tmin <- vector("list", length(topoterra_stack))
for(i in 1:length(topoterra_stack)){
  aet[[i]] <- topoterra_stack[[i]][[1]]
  def[[i]] <- topoterra_stack[[i]][[2]]
  tmax[[i]] <- topoterra_stack[[i]][[3]]
  tmin[[i]] <- topoterra_stack[[i]][[4]]
  
}
#calculate normal
aet_normal <- rast(aet) %>% median()
def_normal<- rast(def) %>% median()
tmax_normal <- rast(tmax) %>% median()
tmin_normal <- rast(tmin) %>% median()

#combine
topoterra_normal <- c(aet_normal, def_normal, tmax_normal, tmin_normal)
names(topoterra_normal) <- names(topoterra_stack[[1]])

#write out
writeRaster(topoterra_normal, "data/merged_output/topoterra_hist_1961-2022.tif",
            overwrite = TRUE, 
            gdal = c("TILED=YES",
                     "COPY_SRC_OVERVIEWS=YES",
                     "COMPRESS=DEFLATE"))

### 2C Normal

#### TerraClimate
terraclimate_files_2c <- list.files("data/cogs", pattern = "2C_[0-9]{4}.tif$", full.names = T)

aet <- terraclimate_files_2c[grep("aet",terraclimate_files_2c)] %>% 
  rast() %>%
  median()
def <- terraclimate_files_2c[grep("def",terraclimate_files_2c)] %>%
  rast() %>%
  median()
tmax <- terraclimate_files_2c[grep("tmax",terraclimate_files_2c)] %>%
  rast() %>%
  median()
tmin <- terraclimate_files_2c[grep("tmin",terraclimate_files_2c)] %>%
  rast() %>%
  median()

terraclimate_normal <- c(aet, def, tmax, tmin) %>%
  project(reference_crs)
names(terraclimate_normal) <- c("aet", "def", "tmax", "tmin")

##  clean terraclim like topoterra
# AET
terraclimate_normal[['aet']][terraclimate_normal[['aet']] < 0] = 0
terraclimate_normal[['aet']][terraclimate_normal[['aet']] > 1000] = 1000
terraclimate_normal[['aet']] <- round(terraclimate_normal[['aet']],0)
# DEF
terraclimate_normal[['def']][terraclimate_normal[['def']] < 0] = 0
terraclimate_normal[['def']][terraclimate_normal[['def']] > 2500] = 2500
terraclimate_normal[['def']] <- round(terraclimate_normal[['def']],0)
# TMAX
terraclimate_normal[['tmax']][terraclimate_normal[['tmax']] < -5] = -5 # This happens to not be needed, but just for completeness
terraclimate_normal[['tmax']][terraclimate_normal[['tmax']] > 35] = 35
terraclimate_normal[['tmax']] <- round(terraclimate_normal[['tmax']], 2)*10^2
# TMIN
terraclimate_normal[['tmin']][terraclimate_normal[['tmin']] < -25] = -25
terraclimate_normal[['tmin']][terraclimate_normal[['tmin']] > 20] = 20
terraclimate_normal[['tmin']] <- round(terraclimate_normal[['tmin']], 2)*10^2

#write out
writeRaster(terraclimate_normal, "data/merged_output/terra_2C_1985-2015.tif",
            overwrite = TRUE, 
            gdal = c("TILED=YES",
                     "COMPRESS=DEFLATE"))

#### TopoTerra
#load topoterra
topoterra_files_2c <- list.files("data/merged_output", pattern = "topoterra_2C", full.names = TRUE)
topoterra_stack <- map(topoterra_files_2c, rast)
#organize into indivudal stacks for each variable
aet <- vector("list", length(topoterra_stack))
def <- vector("list", length(topoterra_stack))
tmax <- vector("list", length(topoterra_stack))
tmin <- vector("list", length(topoterra_stack))
for(i in 1:length(topoterra_stack)){
  aet[[i]] <- topoterra_stack[[i]][[1]]
  def[[i]] <- topoterra_stack[[i]][[2]]
  tmax[[i]] <- topoterra_stack[[i]][[3]]
  tmin[[i]] <- topoterra_stack[[i]][[4]]
  
}
#calculate normal
aet_normal <- rast(aet) %>% median()
def_normal<- rast(def) %>% median()
tmax_normal <- rast(tmax) %>% median()
tmin_normal <- rast(tmin) %>% median()

#combine
topoterra_normal <- c(aet_normal, def_normal, tmax_normal, tmin_normal)
names(topoterra_normal) <- names(topoterra_stack[[1]])

#write out
writeRaster(topoterra_normal, "data/merged_output/topoterra_2C_1985-2015.tif",
            overwrite = TRUE, 
            gdal = c("TILED=YES",
                     "COPY_SRC_OVERVIEWS=YES",
                     "COMPRESS=DEFLATE"))
