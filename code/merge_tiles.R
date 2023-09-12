library(terra)
library(tidyverse)

setwd('C:/Users/PC/Desktop/tyler_working')


list('def','aet','tmin','tmax') %>% 
  walk(\(var){
    grep(var,list.files('GIDs_output_Zoran/gids_output/'), value = T) %>% 
      map(\(tile){
        rast(paste0('GIDs_output_Zoran/gids_output/',tile))
        }) %>% 
      sprc() %>% 
      merge() %>% 
      writeRaster(paste0('data/ds-completed/',var,'_1981-2010.tiff'))
  })
  

