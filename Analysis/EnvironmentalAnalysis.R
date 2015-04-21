library(raster)
library(rgdal)
library(geoR)

# notes from discussion with Amy 12/17
# two things:
#1. Why is latitude a better predictor?
#	a. maybe it's some multivariate thing (see if climate PCs are better than lat or single cliamte var alone)
#	b. maybe need different time averages

#2. Average selective environment
#	

### Occurences
	# Import
	occCard <- read.csv("~/Google Drive/cardinalisENM/data files/all.records.aug.31.csv")
	
	# Convert to SpatialPoints and exclude absences
	occCard <- SpatialPoints(coords = occCard[occCard$PRESABS == 1, c("Longitude", 
		"Latitude")], proj4string = CRS("+proj=longlat +ellps=WGS84"))
	
	# Reproject to match raster layers
	occCard <- spTransform(occCard, CRSobj = CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

### 8 environmental layers used in modeling
	### GET 8 LAYERS FOR FOCAL POPS!
	pred <- stack(paste("~/Google Drive/CardLocalAdaptation/Data/Climate/Raw bios (untransformed)/", c("bio2.grd", "bio3.grd", "bio4.grd", "bio10.grd", "bio11.grd", "bio12.grd", "bio14.grd", "bio15.grd"), sep = ""))

	logPred <- stack(paste("~/Google Drive/CardLocalAdaptation/Data/Climate/Predictor set (log transformed)/", c("bio2.grd", "bio3.grd", "bio4.grd", "bio10.grd", "bio11.grd", "bio12.grd", "bio14.grd", "bio15.grd"), sep = ""))

### Get 8 layers for focal sites
	popenv <- read.csv("~/Google Drive/CardLocalAdaptation/Data/PopulationEnvData.csv")

	# Convert to SpatialPoints
	tmp <- SpatialPoints(coords = popenv[, c("Lon", "Lat")], proj4string = CRS("+proj=longlat +ellps=WGS84"))

	# Reproject to match raster layers
	tmp <- spTransform(tmp, CRSobj = CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

	extract(pred, tmp, method = "bilinear")
	extract(logPred, tmp, method = "bilinear")
### ENM predicted suitability
	
	# Note this currently just from Random Forest model - maybe better to use consensus model or whatever?
	load("~/Google Drive/cardinalisENM/R objects/RF.prob.ave.img")
	# use suitability as weight?
	wts <- extract(RF.prob.ave, occCard)

### PC of 8 climate variables at sites
	climPC <- prcomp(predCard, retx = TRUE, center = TRUE, scale. = TRUE)


###	
### 1. Is latitude complex composite of environmental factors?
###
library(e1071)
occCard <- read.csv("~/Google Drive/cardinalisENM/data files/all.records.aug.31.csv")
absCard <- subset(occCard, occCard$PRESABS == 0)
absCard <- subset(absCard, select = c("Latitude", "Elevation", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15"))

occCard <- subset(occCard, occCard$PRESABS == 1)
occCard <- subset(occCard, select = c("Latitude", "Elevation", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15"))

# Training set
trainset <- sample(1:nrow(occCard), floor(nrow(occCard) / 2))
testset <- c(1:nrow(occCard))[-trainset]

fit <- svm(Latitude ~ Elevation + bio3 + bio4 + bio10 + bio11 + bio12 + bio14 + bio15, data = occCard[trainset, ])

with(occCard[trainset, ], plot(Latitude, predict(fit)))
with(occCard[trainset, ], cor(Latitude, predict(fit)))
with(occCard[testset, ], plot(Latitude, predict(fit, newdata = occCard[testset, -1])))
with(occCard[testset, ], cor(Latitude, predict(fit, newdata = occCard[testset, -1])))

with(absCard, plot(Latitude, predict(fit, newdata = absCard)))
with(absCard, cor(Latitude, predict(fit, newdata = absCard)))

# could test whether this is nonrandom by reperforming with latitudes randomly permuted

occCard$Latitude <- sample(occCard$Latitude, nrow(occCard))
fit <- svm(Latitude ~ Elevation + bio3 + bio4 + bio10 + bio11 + bio12 + bio14 + bio15, data = occCard[trainset, ])

with(occCard[trainset, ], plot(Latitude, predict(fit)))
with(occCard[trainset, ], cor(Latitude, predict(fit)))
with(occCard[testset, ], plot(Latitude, predict(fit, newdata = occCard[testset, -1])))
with(occCard[testset, ], cor(Latitude, predict(fit, newdata = occCard[testset, -1])))


###	
### 2. Does average selective environment change with latitude?
###
	predCard <- extract(pred.dom, occCard)
	plot(occCard$Latitude, predCard[, 1])
	plot(occCard$Latitude, predCard[, 2])
	plot(occCard$Latitude, predCard[, 3])
	plot(occCard$Latitude, predCard[, 4]) # temp of warmest quarter
	plot(occCard$Latitude, predCard[, 5])
	plot(occCard$Latitude, predCard[, 6])
	plot(occCard$Latitude, predCard[, 7])
	plot(occCard$Latitude, predCard[, 8])

### For every presence, calculate nearby spatial variation (coef of var) in suitability
	tmp <- extract(RF.prob.ave, occCard, buffer = 1e4)
	suitCV <- sapply(tmp, function(X) var(X, na.rm = T) / mean(X, na.rm = T))
	plot(occCard$Latitude, suitCV, log = "y")
	fit <- lm(log(suitCV) ~ occCard$Latitude)

	tmp <- extract(pred.dom.wt[[1]], occCard, buffer = 1e4)
	suitBio1 <- sapply(tmp, function(X) var(X, na.rm = T) / mean(X, na.rm = T))
	plot(occCard$Latitude, suitBio1, log = "y")
	fit <- lm(log(suitBio1) ~ occCard$Latitude)

# Variogram or something
	
0.00057485 for buffer = 1e3
	semivarProb <- function(prob, coord, buffer)
	{
		# suitability of occurence points
		probOcc <- extract(prob, coord)

		# suitability of buffer of 1 km
		probBuff <- extract(prob, coord, buffer = buffer)

		# remove NAs
		probBuff <- lapply(probBuff, function(X) X[!is.na(X)])
	
		# remove focal suitability from buffer
		probBuff <- mapply(function(X, Y) X[-which(X == Y)[1]], probBuff, probOcc)

		# calculate semivariance
		out <- mapply(function(X, Y) (1 / (2 * length(X))) * sum((X - Y) ^ 2), probBuff, 
			probOcc)
		
		return(out)
	}

dists <- seq(1e3, 1e5, l = 1e2)
svRF_south <- sapply(dists, function(X) semivarProb(RF.prob.ave, south, buffer = X))
svRF_center <- sapply(dists, function(X) semivarProb(RF.prob.ave, center, buffer = X))
svRF_north <- sapply(dists, function(X) semivarProb(RF.prob.ave, north, buffer = X))
plot(dists, apply(svRF_south, 2, median, na.rm = T), type = "l", ylim = c(0, 0.05))
points(dists, apply(svRF_center, 2, median, na.rm = T), type = "l")
points(dists, apply(svRF_north, 2, median, na.rm = T), type = "l")

for (i in 1:nrow(svRF)) 
plot(dists, svRF[16, ], type = "l", col = "grey", log = "x")
# do by region...
occCard$Latitude



file <- system.file("external/spain.grd", package="usdm")
 
r <- raster(ncols = 100, nrows = 100)
r[] <- rnorm(1e4) 
plot(r) # plot the first RasterLayer in r
 
v1 <- Variogram(r) # compute the sample variogram for the first layer in r
 
v2 <- Variogram(r, cutoff = 30, lag = 2) # specify the lag and cutoff parameters
plot(v2, cloud = T)

RF.prob.ave[1:10, 1:15]