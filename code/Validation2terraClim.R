library(terra)
library(tidyverse)
library(foreach)
library(doParallel)

 ### First Validate Historical data by resampling topoterra to 4km

hist_years <- 1961:2022
hist_residual_list <- foreach(i = hist_years, .packages = c("terra", "tidyverse")) %do% {
  #load topoterra
  topoterra <- rast(paste0("data/merged_output/topoterra_hist_",i,".tif"))
  
  #load terraclim
  terraclim_file <- list.files("data/tile_templates", pattern = paste0("hist_",i,".tif$"), full.names = TRUE)
  terraclim <- rast(terraclim_file)
  ##  clean terraclim like topoterra
  # AET
  terraclim[['aet']][terraclim[['aet']] < 0] = 0
  terraclim[['aet']][terraclim[['aet']] > 1000] = 1000
  terraclim[['aet']] <- round(terraclim[['aet']],0) 
  # DEF
  terraclim[['def']][terraclim[['def']] < 0] = 0
  terraclim[['def']][terraclim[['def']] > 2500] = 2500
  terraclim[['def']] <- round(terraclim[['def']],0)
  # TMAX
  terraclim[['tmax']][terraclim[['tmax']] < -5] = -5 # This happens to not be needed, but just for completeness
  terraclim[['tmax']][terraclim[['tmax']] > 35] = 35
  terraclim[['tmax']] <- round(terraclim[['tmax']], 2)*10^2
  # TMIN
  terraclim[['tmin']][terraclim[['tmin']] < -25] = -25
  terraclim[['tmin']][terraclim[['tmin']] > 20] = 20
  terraclim[['tmin']] <- round(terraclim[['tmin']], 2)*10^2
  
  #resample topoterra to 4km
  topoterra_4km <- resample(topoterra, terraclim, threads = TRUE)
  
  #build topoterra vs terraclim residuals
  residuals <- terraclim - topoterra_4km
  
  #plot histograms of residuals using ggplot, facet by variable name
  df <- residuals %>%
    as.data.frame() %>%
    pivot_longer(cols = everything())
  
  residual_histograms <- df %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 100) +
    facet_wrap(~name, scales = "free") +
    ggtitle(paste0(i, " TopoTerra vs TerraClim Historical Residuals")) +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 14))
  
  final_list <- list(df = df, residual_histograms = residual_histograms)
  }
# save list
save(hist_residual_list, file = "validation/hist_residual_list.RData")
### now perform for 2c

c2_years <- 1985:2015
c2_residual_list <- foreach(i = c2_years, .packages = c("terra", "tidyverse")) %do% {
  #load topoterra
  topoterra <- rast(paste0("data/merged_output/topoterra_2C_",i,".tif"))
  
  #load terraclim
  terraclim_file <- list.files("data/tile_templates", pattern = paste0("2C_",i,".tif$"), full.names = TRUE)
  terraclim <- rast(terraclim_file)
  ##  clean terraclim like topoterra
  # AET
  terraclim[['aet']][terraclim[['aet']] < 0] = 0
  terraclim[['aet']][terraclim[['aet']] > 1000] = 1000
  terraclim[['aet']] <- round(terraclim[['aet']],0) 
  # DEF
  terraclim[['def']][terraclim[['def']] < 0] = 0
  terraclim[['def']][terraclim[['def']] > 2500] = 2500
  terraclim[['def']] <- round(terraclim[['def']],0)
  # TMAX
  terraclim[['tmax']][terraclim[['tmax']] < -5] = -5 # This happens to not be needed, but just for completeness
  terraclim[['tmax']][terraclim[['tmax']] > 35] = 35
  terraclim[['tmax']] <- round(terraclim[['tmax']], 2)*10^2
  # TMIN
  terraclim[['tmin']][terraclim[['tmin']] < -25] = -25
  terraclim[['tmin']][terraclim[['tmin']] > 20] = 20
  terraclim[['tmin']] <- round(terraclim[['tmin']], 2)*10^2
  
  #resample topoterra to 4km
  topoterra_4km <- resample(topoterra, terraclim, threads = TRUE)
  
  #build topoterra vs terraclim residuals
  residuals <- terraclim - topoterra_4km
  
  #plot histograms of residuals using ggplot, facet by variable name
  df <- residuals %>%
    as.data.frame() %>%
    pivot_longer(cols = everything())
  
  residual_histograms <- df %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 100) +
    facet_wrap(~name, scales = "free") +
    ggtitle(paste0(i, " TopoTerra vs TerraClim 2C Residuals")) +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 14))
  
  final_list <- list(df = df, residual_histograms = residual_histograms)
}

#save
save(c2_residual_list, file = "validation/2c_residual_list.RData")

### now perform for 30 year normals- historical ###
#load topoterra normal
topoterra_normal <- rast("data/merged_output/topoterra_hist_1961-2022.tif")
#terraclim_normal 
terraclim_normal <- rast("data/merged_output/terra_hist_1961-2022.tif")


#resample topoterra_normal to terraclim 
topoterra_normal_4k <- resample(topoterra_normal, terraclim_normal, threads = T)

  #build topoterra vs terraclim residuals
  residuals_normal <- terraclim_normal - topoterra_normal_4k
  
  #plot histograms of residuals using ggplot, facet by variable name
  df <- residuals_normal %>%
    as.data.frame() %>%
    pivot_longer(cols = everything())
  
  residual_histograms <- df %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 100) +
    facet_wrap(~name, scales = "free") +
    ggtitle(paste0(i, " TopoTerra vs TerraClim Historical Normals Residuals")) +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 14))
  
  historical_normal_list <- list(df = df, residual_histograms = residual_histograms)
  
  ### repeat for 2C ###
  topoterra_normal <- rast("data/merged_output/topoterra_2C_1985-2015.tif")

#terraclim_normal cleaned to match topoterra

terraclim_normal <- rast("data/merged_output/terra_2C_1985-2015.tif")

#resample topoterra_normal to terraclim
topoterra_normal_4k <- resample(topoterra_normal, terraclim_normal, threads = T)

#build topoterra vs terraclim residuals
residuals_normal <- terraclim_normal - topoterra_normal_4k

#plot histograms of residuals using ggplot, facet by variable name
df <- residuals_normal %>%
  as.data.frame() %>%
  pivot_longer(cols = everything())

residual_histograms <- df %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~name, scales = "free") +
  ggtitle(paste0(i, " TopoTerra vs TerraClim 2C Normals Residuals")) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14))

c2_normal_list <- list(df = df, residual_histograms = residual_histograms)

#combine and save
normal_residual_list <- list(historical_normal_list = historical_normal_list, c2_normal_list = c2_normal_list)
save(normal_residual_list, file = "validation/normal_residual_list.RData")
