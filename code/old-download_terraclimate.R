# 'http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data_plus2C/TerraClimate_2c_vpd_2015.nc'
# download.file(urls, "C:/Users/PC/Desktop/tyler_working/GitHub/topo-downscaling/data/terra_annual/TerraClimate_2c_vpd_2015.nc")

library(purrr)
library(terra)
library(dplyr)
library(tidyr)
library(purrr)

# Climate variables to download
variables <- c('def','aet','tmin','tmax')

# For annualized 
years <- 1985:2015
thredds_2C_base <- paste0('http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data_plus2C/TerraClimate_2c')

years <- 1960:2022
thredds_hist_base <- "http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/data/TerraClimate"

# Build urls based on the base url and variables
urls <- expand.grid(thredds_hist_base, variables, years) |> 
  unite('combo', Var1, Var2, Var3, sep = '_') |> 
  mutate(combo = paste0(combo, '.nc')) 

urls <- split(urls, seq(nrow(urls))) 

urls |> 
  walk(\(url){
    
    # Create a filename to save the netcdf to
    out_file <- paste0('../data/climate_inputs/',sub('.*/data_plus2C/','', url))
    
    # Download and save
    download.file(url$combo, out_file, mode = 'wb')
    
    # Read in as a SpatRast
    nc_as_rast <- rast(out_file)
    
    # Summarize the months to annual
    
    # Extact variable name and year from filename
    variable <- sub('_[0-9]{4}.nc','', sub('.*/TerraClimate_2c_','', url))
    year <- sub('.nc','', sub('.*/TerraClimate_2c_[a-z]{3,4}_','', url))
    
    if (variable %in% c('def','aet')) {
      summ_rast <- sum(nc_as_rast)
    } else {
      summ_rast <- mean(nc_as_rast)
    }
    
    # Save summarized annual as a tiff
    writeRaster(summ_rast, paste0('../data/climate_inputs/',variable,'_terra_2C_',year,'.tif'))
    
    # Delete the netCDF file
    unlink(out_file)
  })

# For 30-yr normals 2C data
thredds_2C_base <- paste0('http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/summaries/TerraClimate2C')

urls <- expand.grid(thredds_2C_base, variables) |> 
  unite('combo', Var1, Var2, sep = '_') |> 
  mutate(combo = paste0(combo, '.nc')) |> 
  split(., seq(nrow(.))) 

urls |> 
  walk(\(url){
    
    # Create a filename to save the netcdf to
    out_file <- paste0('../data/climate_inputs/',sub('.*/TerraClimate2C_','', url))
    
    # Download and save
    download.file(url$combo, out_file, mode = 'wb')
    
    # Read in as a SpatRast
    nc_as_rast <- rast(out_file)
    
    # Summarize the months to annual
    
    # Extact variable name and year from filename
    variable <- sub('*.nc','', sub('.*/TerraClimate2C_','', url))
    
    if (variable %in% c('def','aet')) {
      summ_rast <- sum(nc_as_rast)
    } else {
      summ_rast <- mean(nc_as_rast)
    }
    
    # Save summarized annual as a tiff
    writeRaster(summ_rast, paste0('../data/climate_inputs/',variable,'_terra_2C.tif'))
    
    # Delete the netCDF file
    unlink(out_file)
})

# For 30-yr normals historical data
thredds_2C_base <- paste0('http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/summaries/TerraClimate19611990')

urls <- expand.grid(thredds_2C_base, variables) |> 
  unite('combo', Var1, Var2, sep = '_') |> 
  mutate(combo = paste0(combo, '.nc')) |> 
  split(., seq(nrow(.))) 

urls |> 
  walk(\(url){
    
    # Create a filename to save the netcdf to
    out_file <- paste0('../data/climate_inputs/',sub('.*/TerraClimate19611990_','', url))
    
    # Download and save
    download.file(url$combo, out_file, mode = 'wb')
    
    # Read in as a SpatRast
    nc_as_rast <- rast(out_file)
    
    # Summarize the months to annual
    
    # Extact variable name and year from filename
    variable <- sub('*.nc','', sub('.*/TerraClimate19611990_','', url))
    
    if (variable %in% c('def','aet')) {
      summ_rast <- sum(nc_as_rast)
    } else {
      summ_rast <- mean(nc_as_rast)
    }
    
    # Save summarized annual as a tiff
    writeRaster(summ_rast, paste0('../data/climate_inputs/',variable,'_terra_1961-1990.tif'))
    
    # Delete the netCDF file
    unlink(out_file)
  })


