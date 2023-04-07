library(sf)
library(fasterize)
library(rgeos)
library(raster)
library(parallel)

cwd.250 <- raster('data/topofire/def_topo_1981_2010.tif')
cwd.4km <- raster('data/terraclimate/def_hist_1981_2010.tif')

cwd.250 <- crop(cwd.250, c(-112.5, -111.5, 36, 36.5))
cwd.4km <- crop(cwd.4km, c(-112.5, -111.5, 36, 36.5))

## Make 4km mask for sampling
mask <- projectRaster(cwd.250, cwd.4km, method='bilinear')
plot(mask)

cell.num <- Which(mask < 10000, cells=T)
xy <- as.data.frame(xyFromCell(mask, cell=cell.num))
xy$cell.num <- cell.num
xy$cwd.4km <- round(extract(cwd.4km, xy[, c('x','y')]), 1)
xy$cwd.250 <- round(extract(cwd.250, xy[, c('x','y')], method='bilinear'), 2)

xy.tmp <- xy

cell.num.raster <- rasterFromXYZ(as.matrix(xy[,c('x', 'y', 'cell.num')]))

cell.num.predict <- Which(cwd.250 < 10000, cells=T)
xy.predict <- as.data.frame(xyFromCell(cwd.250, cell=cell.num.predict))
xy.predict$cwd.250 <- round(extract(cwd.250, as.matrix(xy.predict[, c('x','y')]), method='bilinear'), 2)
xy.predict$cell.num <- extract(cell.num.raster, xy.predict[, c('x','y')], method='simple')
xy.predict$cwd.predict <- -99

grca.raster <- crop(mask, c(-113, -111, 35.5, 37))
grca.raster <- extend(grca.raster, mask)
grca.cell.num <- Which(grca.raster >= 0, cells=T)

radius <- 25

##the.function <- function(the.pixel) {
for (i in 1:length(grca.cell.num)) {

	the.pixel <- grca.cell.num[i]
	xy$raster <- 0

	xy[xy$cell.num == the.pixel,]$raster <- 1

	raster.tmp <- rasterFromXYZ(xy[,c('x', 'y', 'raster')])

	raster.tmp[raster.tmp == 0] <- NA

		the.origin <- SpatialPoints(xy.tmp[xy.tmp$cell.num == the.pixel, c('x', 'y')], proj4string=crs(cwd.4km))
		searchLimit <- 6 # the maximum distance in raster units from lakes/boundaries
		polyB <- gBuffer(the.origin, width = searchLimit)
		raster.tmp <- crop(raster.tmp, extent(polyB))

	raster.tmp.dist <- distance(raster.tmp)

	sample.dist <- raster.tmp.dist
	sample.dist[raster.tmp.dist <= radius * 1000] <- 1
	sample.dist[raster.tmp.dist > radius * 1000] <- 0

	cell.num.sample <- Which(sample.dist == 1, cells=T)
	xy.sample <- as.data.frame(xyFromCell(sample.dist, cell.num.sample))
	xy.sample$cwd.250 <- round(extract(cwd.250, xy.sample[, c('x','y')], method='bilinear'), 2)
	xy.sample$cwd.4km <- round(extract(cwd.4km, xy.sample[, c('x','y')]),2)
	xy.sample <- na.omit(xy.sample)

	the.glm <- glm(cwd.4km ~ cwd.250, data=xy.sample, family='quasipoisson')

	xy.predict[which(xy.predict$cell.num == the.pixel),]$cwd.predict <- round(predict(the.glm, newdata=xy.predict[which(xy.predict$cell.num == the.pixel),], type='response'), 0)

##	cwd.predict <- round(predict(the.glm, newdata=xy.predict[which(xy.predict$cell.num == the.pixel),], type='response'), 0)
##	return(cwd.predict)
}

## radius <- 25
## cl <- makeCluster(detectCores()-2)
## 	clusterExport(cl, list('the.function', 'xy', 'xy.tmp', 'cwd.4km', 'cwd.250', 'radius', 'xy.predict', 'cwd.predict'))
##	clusterEvalQ(cl, {
##		library(raster)
##		library(rgeos)})
##	results <- clusterApply(cl, grca.cell.num, the.function)
## stopCluster(cl)


## the.subset <- xy.predict[which(xy.predict$cell.num %in% grca.cell.num),]
## the.subset$cwd.predict <- unlist(results)
the.subset <- subset(xy.predict, cwd.predict != -99)

grca.predict.cwd <- rasterFromXYZ(the.subset[c('x', 'y', 'cwd.predict')])
plot(grca.predict.cwd)

grca.cwd.4km <- mask(cwd.4km, grca.raster)
grca.cwd.4km <- crop(grca.cwd.4km, grca.predict.cwd)

par(mfrow=c(1,2))
plot(grca.cwd.4km)
plot(grca.predict.cwd)

writeRaster(grca.predict.cwd, 'output/def/seans_example.tif')

plot(the.subset$cwd.predict, the.subset$cwd.250)



