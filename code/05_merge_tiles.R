library(terra)
library(tidyverse)


clim_vars <- c('aet','def','tmax','tmin')

# A 'suffix' or naming convention common to all files to be downscaled
coarse_name <- 'terra'

# Time periods
times <- c('1961-1990','2C_1985-2015') #,,paste0('2C_',1985:2015)

out_dir <- 'merged_output'

if (!file.exists(paste0('../data/',out_dir))) {
  dir.create(paste0('../data/',out_dir))
}



as.list(times) %>% 
  walk(\(time){
    clim_vars %>% 
      walk(\(var){
        grep(paste0(var,'_',time),list.files('../data/gids_output/', full.names = T), value = T) %>% 
          map(., rast) %>% 
          sprc() %>% 
          merge() %>% 
          writeRaster(paste0('../data/',out_dir,'/',var,'_',time,'.tif'))
      })
  })

   