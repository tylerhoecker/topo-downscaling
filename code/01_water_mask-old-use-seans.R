# Code to create a "water mask" - a layer which can be used to mask
# climate grids so that values over bodies of water are NA
# This requires the Landfire 2020 BPS data to be on your local machine: https://www.landfire.gov/metadata/lf2020/CONUS/LC20_BPS_220.html

library(terra)

lf_2020 <- rast('data/LC20_BPS_220.tif')

# Aggregate to a coarser resolution, roughly the same as TopoFire (~270 but actually 220 when reprojected to UTM)
# But should be something like: 220/30 = 7.3. This will get it close, and 
# make re-projecting faster, then can resample to match grid and resolution exactly.
lf_2020_coarse <- terra::aggregate(lf_2020, fact = 7, fun = 'modal')

# Skip these steps to do this for entire Landfire coverage 
#-------------------------------------------------------------------------------
# Start by cropping down the LF rast, which is larger than topoclimate data and our study area
topo_area <- rast(rast(file.path('data/climate_inputs/def_topo_hist_1981-2010.tif')))
# Just use the extent to save time
topo_extent <- ext(topo_area)
# Transform to same as Landfire
topo_ext_5070 <- project(topo_extent, from = crs(topo_area), to = crs(lf_2020_coarse))
# Now crop Landire
lf_rast_west <- crop(lf_2020_coarse, topo_ext_5070)
# ------------------------------------------------------------------------------

# Reclassify landfire
lf_rast_bin <- classify(lf_rast_west, cbind(c(-9999,11),NA))
lf_rast_bin <- classify(lf_rast_bin, cbind(c(0:9999),1))

# ------------------------------------------------------------------------------
lf_water_mask <- project(lf_rast_bin, "EPSG:32611")
writeRaster(lf_water_mask, 'lf_water_mask_32611.tiff', overwrite = T)

# Or do re-projection to UTM Zone 11 EPSG:3741 in QGIS because its WAY faster!
writeRaster(lf_rast_bin, 'lf_water_mask_5070.tiff', overwrite = T)
# ------------------------------------------------------------------------------




# Write out


