library(terra)
library(tidyverse)

# Time periods
times <- paste0('2C_',1985:2015) #,,'1961-1990',paste0('2C_',1985:2015)

out_dir <- 'merged_output'

if (!file.exists(paste0('../data/',out_dir))) {
  dir.create(paste0('../data/',out_dir))
}

as.list(times) %>% 
  walk(\(time){
    grep(time,list.files('../data/gids_output', full.names = T), value = T) %>% 
      map(., rast) %>% 
      sprc() %>% 
      merge() %>% 
      writeRaster(paste0('../data/',out_dir,'/','topoterra_',time,'.tif'))
  })

   