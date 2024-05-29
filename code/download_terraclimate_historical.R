# 'http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data_plus2C/TerraClimate_2c_vpd_2015.nc'
# download.file(urls, "C:/Users/PC/Desktop/tyler_working/GitHub/topo-downscaling/data/terra_annual/TerraClimate_2c_vpd_2015.nc")

library(purrr)
library(terra)
library(dplyr)
library(tidyr)
library(purrr)

# Climate variables to download
variables <- c("tmax","tmin") #c('def','aet','tmin','tmax')

# For annualized 
years <- 1985:2015
thredds_base <- paste0('http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data_plus2C/TerraClimate_2c')

#years <- 1960:2022
#thredds_base <- "http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/TerraClimate"

# Build urls based on the base url and variables
urls <- expand.grid(thredds_base, variables, years) |> 
  unite('combo', Var1, Var2, Var3, sep = '_') |> 
  mutate(combo = paste0(combo, '.nc')) 

urls <- split(urls, seq(nrow(urls))) 

urls |> 
  walk(\(url){
    
    # Create a filename to save the netcdf to
    out_file <- paste0('data/climate_inputs/',sub('.*/data.*/','', url))
    
    # Download and save
    download.file(url$combo, out_file, mode = 'wb')
    
    # Read in as a SpatRast
    nc_as_rast <- rast(out_file)
        
    # Extact variable name and year from filename
    variable <- varnames(nc_as_rast)
    year <- sub('.nc','', sub('.*/TerraClimate_.+[a-z]{3,4}_','', url))
    period <- ifelse(
      grepl("2c", url),
      "2C",
      "hist"
    )
    
    # Summarize monthly data to appropriate annual statistic
    if (variable %in% c('def','aet')) {
      summ_rast <- sum(nc_as_rast)
    } 
    
    if (variable == "tmax") {
      summ_rast <- nc_as_rast[[7]]
    } 
    
    if (variable == "tmin") {
      summ_rast <- nc_as_rast[[1]]
    }

    # Save summarized annual as a tiff
    writeRaster(summ_rast, paste0('data/climate_inputs/',variable,'_terra_',period,'_',year,'.tif'))
    
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


