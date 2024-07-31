# This script validates topoTerra data against station data. 

# The script reads in the station data and the topoTerra data,
# then compares the two datasets.
library(tidyverse)
library(terra)
library(data.table)
library(sf)
# Load the data
stations <- vect("validation/station_data/complete_station_data.gpkg")

# clean stations like topoterra
# TMAX
stations$TMAX[stations$TMAX < -5] = -5 # This happens to not be needed, but just for completeness
stations$TMAX[stations$TMAX > 35] = 35
stations$TMAX <- round(stations$TMAX, 2)*10^2
# TMIN
stations$TMIN[stations$TMIN < -25] = -25
stations$TMIN[stations$TMIN > 20] = 20
stations$TMIN <- round(stations$TMIN, 2)*10^2

years <- 1961:2022
topoterra <- list.files("data/merged_output/",
                        pattern = "hist_[0-9]{4}.tif",
                        full.names = TRUE) %>%
  map(terra::rast)

# Extract the raster data to the station locations by year
yearly_station_data <- vector("list", length(years))
for( year_i in seq_along(years)){
    year <- years[year_i]

    yearly_station_data[[year_i]] <- stations %>%
        subset(YEAR == year, NSE = TRUE) %>%
        terra::extract(topoterra[[year_i]], .) %>%
        cbind(stations %>%
                  subset(YEAR == year, NSE = TRUE, select = c(station_id, YEAR))) %>%
        terra::subset(!is.na(station_id), select = c(station_id, YEAR, tmax, tmin), NSE =TRUE)
}

yearly_station_data <- rbindlist(yearly_station_data)

# join the yearly station_data to the original station data
validation_stations <- stations %>%
  as.data.frame(geom = "WKT") %>%  
  left_join(yearly_station_data, by = c("station_id", "YEAR")) %>%
  rename(
    tmax_topoterra = "tmax",
    tmin_topoterra = "tmin",
    tmax_station = "TMAX",
    tmin_station = "TMIN"
    
  )

# save the data
saveRDS(validation_stations, "validation/station_data/validation_stations.rds")

# now we can use the validation_stations data to compare the station data to the topoterra data

# The script reads in the station data and the topoTerra data,
# then compares the two datasets.

validation_stations <- readRDS("validation/station_data/validation_stations.rds")

# Calculate the difference between the station data and the topoterra data
validation_stations <- validation_stations %>%
  mutate(
    tmax_diff = tmax_station - tmax_topoterra,
    tmin_diff = tmin_station - tmin_topoterra
  )

  # plot observed vs predicted(topoterra)
observed_predicted_plot_tmax <- validation_stations %>%
    ggplot(aes(x = tmax_station, y = tmax_topoterra)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
        x = "Observed (Station) TMAX",
        y = "Predicted (TopoTerra) TMAX",
        title = "Observed vs Predicted TMAX"
    ) +
    theme_minimal()
ggsave( "validation/figures/observed_predicted_plot_tmax.jpg", observed_predicted_plot_tmax)

observed_predicted_plot_tmin <- validation_stations %>%
    ggplot(aes(x = tmin_station, y = tmin_topoterra)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
        x = "Observed (Station) TMIN",
        y = "Predicted (TopoTerra) TMIN",
        title = "Observed vs Predicted TMIN"
    ) +
    theme_minimal()
ggsave( "validation/figures/observed_predicted_plot_tmin.jpg", observed_predicted_plot_tmin)

#calculate the correlation between the station data and the topoterra data for each station
correlation <- validation_stations %>%
  group_by(station_id) %>%
  summarise(
    tmax_correlation = cor(tmax_station, tmax_topoterra, use = "complete.obs"),
    tmin_correlation = cor(tmin_station, tmin_topoterra, use = "complete.obs")
  ) %>%
  # if NA, set to 1
    mutate(
        tmax_correlation = ifelse(is.na(tmax_correlation), 1, tmax_correlation),
        tmin_correlation = ifelse(is.na(tmin_correlation), 1, tmin_correlation)
    )

# bind to stations data with coordinates
correlation_sf <- stations %>%
  as.data.frame(geom = "WKT") %>%
  select(station_id, geometry) %>%
  distinct() %>%
  left_join(correlation, by = c("station_id" = "station_id")) %>%
  vect(., geom = c("geometry"), crs = crs(stations)) %>%
  st_as_sf()

#plot the correlation on a map
#load western states
states <- st_read("data/western_states/western_states.shp") %>%
  st_transform(crs = crs(correlation_sf))
correlation_sf <- st_intersection(correlation_sf, states)
correlation_plot_tmax <-  ggplot() +
  geom_sf( fill = "grey70", color = "black", data = states) +
  geom_sf(aes(color = tmax_correlation), data = correlation_sf) +
  scale_color_viridis_b(option = "magma") +
  labs(
    title = "TMAX Correlation",
    color = "Correlation"
  ) +
  theme_minimal()
ggsave( "validation/figures/correlation_plot_tmax.jpg", correlation_plot_tmax)

correlation_plot_tmin <-  ggplot() +
    geom_sf( fill = "grey70", color = "black", data = states) +
    geom_sf(aes(color = tmin_correlation), data = correlation_sf) +
    scale_color_viridis_b(option = "magma") +
    labs(
        title = "TMIN Correlation",
        color = "Correlation"
    ) +
    theme_minimal()
ggsave( "validation/figures/correlation_plot_tmin.jpg", correlation_plot_tmin)
