library(raster)
library(rgdal)

### Occurences
	# Import
	occCard <- read.csv("~/Google Drive/cardinalisENM/data files/herb.dat.lnpreds.csv")
	
	# Convert to SpatialPoints and exclude absences
	occCard <- SpatialPoints(coords = occCard[occCard$PRESABS == 1, c("Longitude", 
		"Latitude")], proj4string = CRS("+proj=longlat +ellps=WGS84"))
	
	# Reproject to match raster layers
	occCard <- spTransform(occCard, CRSobj = CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
	
### 8 environmental layers used in modeling
	
	# This doesn't work because file path is wrong
	# load("/Users/cdmuir/Google Drive/cardinalisENM/R objects/wc0.5/pred.dom.Rda")
	
	# This recreates the same thing with correct local path
	pred.dom <- stack(paste("~/Google Drive/cardinalisENM/R objects/wc0.5/", c("bio2.grd", 
		"bio3.grd", "bio4.grd", "bio10.grd", "bio11.grd", "bio12.grd", "bio14.grd", 
		"bio15.grd"), sep = ""))

### ENM predicted suitability

	load("~/Google Drive/cardinalisENM/R objects/RF.prob.ave.img")

### Overlay to get predictors weighted by suitability
	
	pred.dom.wt <- overlay(pred.dom, RF.prob.ave, fun = function(x, y) {(x * y)})

### For every presence, calculate nearby spatial variation (coef of var) in suitability
	tmp <- extract(RF.prob.ave, occCard, buffer = 1e4)
	suitCV <- sapply(tmp, function(X) var(X, na.rm = T) / mean(X, na.rm = T))
	plot(occCard$Latitude, suitCV, log = "y")
	fit <- lm(log(suitCV) ~ occCard$Latitude)

	tmp <- extract(pred.dom.wt[[1]], occCard, buffer = 1e4)
	suitBio1 <- sapply(tmp, function(X) var(X, na.rm = T) / mean(X, na.rm = T))
	plot(occCard$Latitude, suitBio1, log = "y")
	fit <- lm(log(suitBio1) ~ occCard$Latitude)
	