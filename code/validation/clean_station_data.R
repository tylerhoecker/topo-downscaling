#script that takes the raw station data and cleans it up for use in the analysis

library(tidyverse)
library(terra)
library(VFS)
library(furrr)
library(data.table)
# Load the data
stations <- read_table("validation/station_data/ghcnd-stations.txt", col_names = FALSE) %>%
  rename(
    station_id = X1,
    latitude = X2,
    longitude = X3,
    elevation = X4,
    state = X5,
    name = X6,
    gsn_flag = X7,
    hcn_crn_flag = X8,
    wmo_id = X9
  )
stations_v <- vect(stations, geom = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")

#load template
template <- rast("data/merged_output/topoterra_hist_1961-2022.tif")[[1]]

# match stations_v to template
stations_v <- project(stations_v, template) %>%
  crop(template)

# extract the station data
stations_US <- extract(template, stations_v, bind = T)

#remove na aet
stations_US <- stations_US[!is.na(stations_US$aet),]

#extract station IDS
station_ids <- stations_US$station_id

#list station data files, then filter by station ids
stations <- list.files("validation/station_data/ghcnd_all/", pattern = ".dly", full.names = T) %>%
  as_tibble() %>%
  filter(str_sub(basename(value), 1, 11) %in% station_ids) %>%
  pull()


#reduce and clean station data
year_range <- 1961:2022
plan(multisession, workers = 20)
valid_stations <- future_map_lgl(stations, \(station){
  # read in the data
  station <- read.dly(station)
  
  filtered_station <- station %>%
    # filter to years
    filter(YEAR %in% year_range)
  
  # check if there are 30 complete years within filtered stations, return NULL if not
  complete_years <- filtered_station %>%
    drop_na(TMIN.VALUE, TMAX.VALUE) %>%
    group_by(YEAR) %>%
    # check if there are 365 or 366 days in the year
    summarise(n = n()) %>%
    filter(n == 365 | n == 366) %>%
    nrow()
  
  if (complete_years < 30) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
plan(sequential)

valid_station_files <- stations[valid_stations]

# create monthly data for each station, then summarize to yearly
station_data <- future_map(valid_station_files, \(station){
  station_name <- basename(station) %>%
    str_remove(".dly")
  # read in the data
  station <- read.dly(station) %>%
  select(all_of(c("YEAR", "MONTH", "TMAX.VALUE", "TMIN.VALUE")))
  
  # filter to years
  filtered_station <- station %>%
    filter(YEAR %in% year_range)

    # filter to complete years
  filtered_station <- filtered_station %>%
    drop_na(TMIN.VALUE, TMAX.VALUE) %>%
    group_by(YEAR) %>%
    filter(n() == 365 | n() == 366)
  
  # create monthly summaries
  monthly_data <- filtered_station %>%
    group_by(YEAR, MONTH) %>%
    summarise(
      TMAX = mean(TMAX.VALUE, na.rm = TRUE),
      TMIN = mean(TMIN.VALUE, na.rm = TRUE)
    )
  
  # summarize to yearly 
  # tmax is tmax during July
  # tmin is tmin during January
  yearly_data <- monthly_data %>%
    group_by(YEAR) %>%
    summarise(
      TMAX = TMAX[MONTH == 7],
      TMIN = TMIN[MONTH == 1]
    )
    # add station_id
    yearly_data$station_id <- station_name
  return(yearly_data)
})

# bind all the station data together
station_data <- rbindlist(station_data)

# join station_data with  stations_US
  yearly_station_data <- as.data.frame(stations_US, geom = "WKT") %>%
  right_join(station_data, by = c("station_id" = "station_id")) %>%
  vect(., geom = c("geometry"), crs = crs(stations_US)) %>%
  terra::subset(subset = !is.na(station_id), select = c(station_id, YEAR, TMAX, TMIN), NSE = TRUE)

# save the data
writeVector(yearly_station_data, "validation/station_data/complete_station_data.gpkg", overwrite = TRUE)
 
      
