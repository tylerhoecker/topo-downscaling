library(sf)
library(fasterize)
library(rgeos)
library(raster)
library(parallel)

cwd.250 <- raster('data/topofire/def_1981_2010.tif')
cwd.4km <- raster('data/terraclimate/def_hist_1981_2010.tif')

cwd.250 <- crop(cwd.250, c(-112, -100, 32, 44))
cwd.4km <- crop(cwd.4km, c(-114, -98, 30, 46))

## Make 4km mask for sampling
mask <- projectRaster(cwd.250, cwd.4km, method='bilinear')
plot(mask)

cell.num <- Which(mask < 10000, cells=T)
xy <- as.data.frame(xyFromCell(mask, cell=cell.num))
xy$cell.num <- cell.num
xy$cwd.4km <- round(extract(cwd.4km, xy[, c('x','y')]), 1)
xy$cwd.250 <- round(extract(cwd.250, xy[, c('x','y')], method='bilinear'), 2)

xy.tmp <- xy

## Get random sample  of 1000 250m pixels to test various radii
set.seed(666)
xy.rand <- xy[sample(1:nrow(xy), 1000, replace=F),]
points(xy.rand$x, xy.rand$y)

the.function <- function(the.row) {
	xy.tmp$raster <- 0

	cell.id <- row.names(xy.rand[the.row,])
	xy.tmp[cell.id,]$raster <- 1
	
	raster.tmp <- rasterFromXYZ(xy.tmp[c('x', 'y', 'raster')])
	raster.tmp[raster.tmp == 0] <- NA

	raster.tmp.dist <- distance(raster.tmp)

	sample.dist <- raster.tmp.dist
	sample.dist[raster.tmp.dist <= radius * 1000] <- 1
	sample.dist[raster.tmp.dist > radius * 1000] <- 0

	xy.tmp$sample <- extract(sample.dist, xy.tmp[,c('x', 'y')])
	xy.sample <- subset(xy.tmp, sample == 1)

	## I don't know if 'quasipoisson' is the correct family, but I do no that there are no predictions
	## less than zero which could occur with a standard Gaussian distribution. Plus it seems to work
	the.glm <- glm(cwd.4km ~ cwd.250, data=xy.sample, family='quasipoisson')

	the.prediction <- round(predict(the.glm, newdata=xy.tmp[cell.id,], type='response'))
	
	return(the.prediction)
}

radius <- 25
cl <- makeCluster(detectCores()-2)

clusterExport(cl, list('the.function', 'xy.tmp', 'xy.rand', 'radius'))
clusterEvalQ(cl, {library(raster)})
results <- clusterApply(cl, 1:nrow(xy.rand), the.function)

stopCluster(cl)

xy.rand$predict.cwd <- unlist(results)
xy.rand.1 <- subset(xy.rand, predict.cwd != 'NA')
xy.rand.1

## Test results vs. 4km TerraClimate
## I'm not sure highest R2 is really that important, but at least this give some indication that the procedure is working

## 25 km (115 pixels in model)
(r2 <- round(summary(lm(cwd.4km ~ predict.cwd, data=xy.rand.1))$r.squared, 2)) ## 0.99
(slope <- round(summary(lm(cwd.4km ~ predict.cwd, data=xy.rand.1))$coef[2,1], 2)) ## 1.00
plot(xy.rand.1$predict.cwd, xy.rand.1$cwd.4km)

## 50 km (451 pixels in model)
(r2 <- round(summary(lm(cwd.4km ~ predict.cwd, data=xy.rand.1))$r.squared, 2)) ## 0.98
(slope <- round(summary(lm(cwd.4km ~ predict.cwd, data=xy.rand.1))$coef[2,1], 2)) ## 1.00
plot(xy.rand.1$predict.cwd, xy.rand.1$cwd.4km)

## 100 km (1803 pixels in model)
(r2 <- round(summary(lm(cwd.4km ~ predict.cwd, data=xy.rand.1))$r.squared, 2)) ## 0.97
(slope <- round(summary(lm(cwd.4km ~ predict.cwd, data=xy.rand.1))$coef[2,1], 2)) ## 1.00
plot(xy.rand.1$predict.cwd, xy.rand.1$cwd.4km)


