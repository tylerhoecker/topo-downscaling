#download GHVCv3 data from the website
#
library(archive)
library(tidyverse)
library(terra)
library(sf)
# download dly files
options(timeout=60*10)
download.file("https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd_all.tar.gz", "validation/station_data/ghcnd_all.tar.gz")
archive_extract("validation/station_data/ghcnd_all.tar.gz", dir = "validation/station_data")

# download metadata
download.file("https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt", "validation/station_data/ghcnd-stations.txt")
