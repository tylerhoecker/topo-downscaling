# This script validates topoTerra data against station data. 

# The script reads in the station data and the topoTerra data,
# then compares the two datasets.
library(tidyverse)
library(terra)
library(data.table)
library(sf)
library(ggpubr)
# Load the data
stations <- vect("validation/station_data/complete_station_data.gpkg")

# clean stations like topoterra
# TMAX
stations$TMAX[stations$TMAX < -5] = -5 # This happens to not be needed, but just for completeness
stations$TMAX[stations$TMAX > 50] = 50
stations$TMAX <- round(stations$TMAX, 2)
# TMIN
stations$TMIN[stations$TMIN < -25] = -25
stations$TMIN[stations$TMIN > 20] = 20
stations$TMIN <- round(stations$TMIN, 2)

years <- 1961:2022
topoterra <- list.files("data/merged_output/",
                        pattern = "hist_[0-9]{4}.tif",
                        full.names = TRUE) %>%
  map(terra::rast)

terraclimate <- list.files("data/tile_templates/",
                           pattern = "hist_[0-9]{4}.tif",
                           full.names = TRUE) %>%
  map(terra::rast)

# Extract the raster data to the station locations by year
yearly_station_data_topoterra <- vector("list", length(years))
yearly_station_data_terraclimate <- vector("list", length(years))
for( year_i in seq_along(years)){
    year <- years[year_i]

    yearly_station_data_topoterra[[year_i]] <- stations %>%
        terra::subset(subset = YEAR == year, NSE = TRUE) %>%
        terra::extract(topoterra[[year_i]], ., method = "bilinear") %>%
        cbind(stations %>%
                 terra::subset(YEAR == year, NSE = TRUE, select = c(station_id, YEAR))) %>%
        base::subset(!is.na(station_id), select = c(station_id, YEAR, tmax, tmin))
    
    yearly_station_data_terraclimate[[year_i]] <- stations %>%
        terra::subset(YEAR == year, NSE = TRUE) %>%
        terra::extract(terraclimate[[year_i]], ., method = "bilinear") %>%
        cbind(stations %>%
                  subset(YEAR == year, NSE = TRUE, select = c(station_id, YEAR))) %>%
        base::subset(!is.na(station_id), select = c(station_id, YEAR, tmax, tmin))
}

yearly_station_data_topoterra <- rbindlist(yearly_station_data_topoterra)

yearly_station_data_terraclimate <- rbindlist(yearly_station_data_terraclimate)

# join the yearly station_data to the original station data
validation_stations_topoterra <- stations %>%
  as.data.frame(geom = "WKT") %>%  
  left_join(yearly_station_data_topoterra, by = c("station_id", "YEAR")) %>%
  rename(
    tmax_topoterra = "tmax",
    tmin_topoterra = "tmin",
    tmax_station = "TMAX",
    tmin_station = "TMIN"
    
  )

  ## terraclimate
validation_stations_terraclimate <- stations %>%
  as.data.frame(geom = "WKT") %>%  
  left_join(yearly_station_data_terraclimate, by = c("station_id", "YEAR")) %>%
  rename(
    tmax_terraclimate = "tmax",
    tmin_terraclimate = "tmin",
    tmax_station = "TMAX",
    tmin_station = "TMIN"
    
  )

# save the data
saveRDS(validation_stations_topoterra, "validation/station_data/validation_stations_topoterra.rds")
saveRDS(validation_stations_terraclimate, "validation/station_data/validation_stations_terraclimate.rds")

# now we can use the validation_stations_topoterra data to compare the station data to the topoterra data

# The script reads in the station data and the topoTerra data,
# then compares the two datasets.

validation_stations_topoterra <- readRDS("validation/station_data/validation_stations_topoterra.rds")
validation_stations_terraclimate <- readRDS("validation/station_data/validation_stations_terraclimate.rds")

# Calculate the difference between the station data and the topoterra data
validation_stations_topoterra <- validation_stations_topoterra %>%
  mutate(
    tmax_station = tmax_station,
    tmin_station = tmin_station,
    tmax_topoterra = tmax_topoterra,
    tmin_topoterra = tmin_topoterra,
    tmax_diff = (tmax_station - tmax_topoterra),
    tmin_diff = (tmin_station - tmin_topoterra)
  )

validation_stations_terraclimate <- validation_stations_terraclimate %>%
  mutate(
    tmax_station = tmax_station,
    tmin_station = tmin_station,
    tmax_diff = (tmax_station - tmax_terraclimate),
    tmin_diff = (tmin_station - tmin_terraclimate)
  )

  # plot observed vs predicted(topoterra)
  pearson_tmax_topoterra <- cor(validation_stations_topoterra$tmax_station, validation_stations_topoterra$tmax_topoterra)
  r2_tmax_topoterra <- summary(lm(tmax_topoterra ~ tmax_station, data = validation_stations_topoterra))$adj.r.squared
observed_predicted_plot_tmax_topoterra <- validation_stations_topoterra %>%
    ggplot(aes(x = tmax_station, y = tmax_topoterra)) +
    geom_point() +
    stat_regline_equation(
        aes(label = paste(after_stat(adj.rr.label), sep = "~~~~")),
        formula = y ~ x,
        label.x = 0.1,
        label.y = 40,
        size = 8
    ) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
        x = "Observed (Station) TMAX (\u00B0C)",
        y = "Predicted (TopoTerra) TMAX (\u00B0C)",
        title = "Observed vs Predicted TMAX (TopoTerra)"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 16))
ggsave( "validation/figures/observed_predicted_plot_tmax_topoterra.jpg", observed_predicted_plot_tmax_topoterra)

pearson_tmax_terra <- cor(validation_stations_terraclimate$tmax_station, validation_stations_terraclimate$tmax_terraclimate)
r2_tmax_terra <- summary(lm(tmax_terraclimate ~ tmax_station, data = validation_stations_terraclimate))$adj.r.squared
observed_predicted_plot_tmax_terraclimate <- validation_stations_terraclimate %>%
    ggplot(aes(x = tmax_station, y = tmax_terraclimate)) +
    geom_point() +
    stat_regline_equation(
        aes(label = paste(after_stat(adj.rr.label), sep = "~~~~")),
        formula = y ~ x,
        label.x = 0.1,
        label.y = 40,
        size = 8
    ) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
        x = "Observed (Station) TMAX (\u00B0C)",
        y = "Predicted (TerraClimate) TMAX (\u00B0C)",
        title = "Observed vs Predicted TMAX (TerraClimate)"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 16))

ggsave( "validation/figures/observed_predicted_plot_tmax_terraclimate.jpg", observed_predicted_plot_tmax_terraclimate)


#plot observed vs predicted(tmin)
pearson_tmin_topoterra <- cor(validation_stations_topoterra$tmin_station, validation_stations_topoterra$tmin_topoterra)
r2_tmin_topoterra <- summary(lm(tmin_topoterra ~ tmin_station, data = validation_stations_topoterra))$adj.r.squared
observed_predicted_plot_tmin_topoterra <- validation_stations_topoterra %>%
    ggplot(aes(x = tmin_station, y = tmin_topoterra)) +
    geom_point() +
    stat_regline_equation(
        aes(label = paste(after_stat(adj.rr.label), sep = "~~~~")),
        formula = y ~ x,
        label.x = -20,
        label.y = 40,
        size = 8
    ) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
        x = "Observed (Station) TMIN (\u00B0C)",
        y = "Predicted (TopoTerra) TMIN (\u00B0C)",
        title = "Observed vs Predicted TMIN (TopoTerra)"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 16))
ggsave( "validation/figures/observed_predicted_plot_tmin_topoterra.jpg", observed_predicted_plot_tmin_topoterra)

pearson_tmin_terra <- cor(validation_stations_terraclimate$tmin_station, validation_stations_terraclimate$tmin_terraclimate)
r2_tmin_terra <- summary(lm(tmin_terraclimate ~ tmin_station, data = validation_stations_terraclimate))$adj.r.squared
observed_predicted_plot_tmin_terraclimate <- validation_stations_terraclimate %>%
    ggplot(aes(x = tmin_station, y = tmin_terraclimate)) +
    geom_point() +
    stat_regline_equation(
        aes(label = paste(after_stat(adj.rr.label), sep = "~~~~")),
        formula = y ~ x,
        label.x = -20,
        label.y = 40,
        size = 8
    ) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
        x = "Observed (Station) TMIN (\u00B0C)",
        y = "Predicted (TerraClimate) TMIN (\u00B0C)",
        title = "Observed vs Predicted TMIN (TerraClimate)"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 16))
ggsave( "validation/figures/observed_predicted_plot_tmin_terraclimate.jpg", observed_predicted_plot_tmin_terraclimate)

# save each r and r^2 in a dataframe with data set and variable
r2_data <- data.frame(
    data_set = c("TopoTerra", "TopoTerra", "TerraClimate", "TerraClimate"),
    variable = c("TMAX", "TMIN", "TMAX", "TMIN"),
    r = c(pearson_tmax_topoterra, pearson_tmin_topoterra, pearson_tmax_terra, pearson_tmin_terra),
    r2 = c(r2_tmax_topoterra, r2_tmin_topoterra, r2_tmax_terra, r2_tmin_terra)
)
write.csv(r2_data, "validation/figures/r2_data.csv", row.names = FALSE)
#calculate the correlation_topoterra between the station data and the topoterra data for each station
correlation_topoterra <- validation_stations_topoterra %>%
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

correlation_terraclimate <- validation_stations_terraclimate %>%
  group_by(station_id) %>%
  summarise(
    tmax_correlation = cor(tmax_station, tmax_terraclimate, use = "complete.obs"),
    tmin_correlation = cor(tmin_station, tmin_terraclimate, use = "complete.obs")
  ) %>%
  # if NA, set to 1
    mutate(
        tmax_correlation = ifelse(is.na(tmax_correlation), 1, tmax_correlation),
        tmin_correlation = ifelse(is.na(tmin_correlation), 1, tmin_correlation)
    )

# bind to stations data with coordinates
correlation_sf_topoterra <- stations %>%
  as.data.frame(geom = "WKT") %>%
  select(station_id, geometry) %>%
  distinct() %>%
  left_join(correlation_topoterra, by = c("station_id" = "station_id")) %>%
  vect(., geom = c("geometry"), crs = crs(stations)) %>%
  st_as_sf()

correlation_sf_terraclimate <- stations %>%
  as.data.frame(geom = "WKT") %>%
  select(station_id, geometry) %>%
  distinct() %>%
  left_join(correlation_terraclimate, by = c("station_id" = "station_id")) %>%
  vect(., geom = c("geometry"), crs = crs(stations)) %>%
  st_as_sf()
#plot the correlation_topoterra on a map
#load western states
states <- st_read("data/western_states/western_states.shp") %>%
  st_transform(crs = crs(correlation_sf_topoterra))
correlation_sf_topoterra <- st_intersection(correlation_sf_topoterra, states)
correlation_sf_terraclimate <- st_intersection(correlation_sf_terraclimate, states)
correlation_plot_tmax_topoterra <-  ggplot() +
  geom_sf( fill = "grey70", color = "black", data = states) +
  geom_sf(aes(color = tmax_correlation), data = correlation_sf_topoterra) +
  scale_color_viridis_b(option = "magma") +
  labs(
    title = "TMAX Correlation (Station vs TopoTerra)",
    color = "Correlation"
  ) +
  theme_minimal() +
    theme(text = element_text(size = 16))
ggsave( "validation/figures/correlation_plot_tmax_topoterra.jpg", correlation_plot_tmax_topoterra)

correlation_plot_tmax_terraclimate <-  ggplot() +
  geom_sf( fill = "grey70", color = "black", data = states) +
  geom_sf(aes(color = tmax_correlation), data = correlation_sf_terraclimate) +
  scale_color_viridis_b(option = "magma") +
  labs(
    title = "TMAX Correlation (Station vs TerraClimate)",
    color = "Correlation"
  ) +
  theme_minimal() +
    theme(text = element_text(size = 16))
ggsave( "validation/figures/correlation_plot_tmax_terraclimate.jpg", correlation_plot_tmax_terraclimate)

correlation_plot_tmin_topoterra <-  ggplot() +
    geom_sf( fill = "grey70", color = "black", data = states) +
    geom_sf(aes(color = tmin_correlation), data = correlation_sf_topoterra) +
    scale_color_viridis_b(option = "magma") +
    labs(
        title = "TMIN Correlation (Station vs TopoTerra)",
        color = "Correlation"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 16))
ggsave( "validation/figures/correlation_plot_tmin_topoterra.jpg", correlation_plot_tmin_topoterra)

correlation_plot_tmin_terraclimate <-  ggplot() +
    geom_sf( fill = "grey70", color = "black", data = states) +
    geom_sf(aes(color = tmin_correlation), data = correlation_sf_terraclimate) +
    scale_color_viridis_b(option = "magma") +
    labs(
        title = "TMIN Correlation (Station vs TerraClimate)",
        color = "Correlation"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 16))
ggsave( "validation/figures/correlation_plot_tmin_terraclimate.jpg", correlation_plot_tmin_terraclimate)

