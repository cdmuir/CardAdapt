
# Libraries

	library(raster)
	library(rgeos)

# Color scheme

	colNorth <- "deepskyblue"
	colSouth <- rgb(255, 131, 0, maxColorValue = 255)

	# Get maps of California and Oregon
  usa <- getData('GADM', country = "USA", level = 1, path = pathDat)
  states <- gSimplify(usa[usa@data$NAME_1 %in% c("California", "Oregon", "Nevada", "Idaho"), ], 
										tol = 0.01, topologyPreserve = TRUE)

  # Load thinned occurences used Climate-Lat analyses and convert to Spatial class
  occ <- read.csv(paste(pathDat, "/thinned1OccLatLon.csv", sep = ""), row.names = NULL)
  coordinates(occ) <- ~Longitude+Latitude
  projection(occ) <- CRS('+proj=longlat');

  # Load focal populations and convert to Spatial class
  load(paste(pathDat, "/focClim.RData", sep = ""))
  foc <- focClim[, c("Longitude", "Latitude")]
  coordinates(foc) <- ~Longitude+Latitude
  projection(foc) <- CRS('+proj=longlat');

  # Clip and map
  clipExtent <- as(extent(floor(min(occ$Longitude) * 10) / 10, 
	   											ceiling(max(occ$Longitude) * 10) / 10, 
	  											floor(min(occ$Latitude) * 10) / 10, 
		  										ceiling(max(occ$Latitude) * 10) / 10), 
			  					 "SpatialPolygons")
  proj4string(clipExtent) <- CRS(proj4string(states))
  states <- gIntersection(states, clipExtent, byid = TRUE)

  # Transform to Albert's Equal Area Projection
  projAEA <- sprintf("+proj=aea +lat_1=%s +lat_2=%s +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
	  								 extent(clipExtent)[3], extent(clipExtent)[4], 
		  							 extent(clipExtent)[1], extent(clipExtent)[2])
  states <- spTransform(states, CRS(projAEA))
  occAEA <- spTransform(occ, CRS(projAEA))
  focAEA <- spTransform(foc, CRS(projAEA))

  # PLOT
  pdf(paste(pathFig, "/Figure_VarSelectBase.pdf", sep = ""), w = 6.5, h = 7, 
	  	useDingbat = F)
  mat <- matrix(c(1, 4, 7,
	  							1, 4, 7,
		  						1, 4, 7,
			  					2, 5, 8,
				  				2, 5, 8,
					  			3, 6, 9,
						  		3, 6, 9), ncol = 3, byrow = T)
  layout(mat)

  # A. Example Climate-Latitude variable

	  ## Panel 1: Map
	  par(mai = rep(0.5, 4))
	  plot(states, main = "A. Climate-Latitude variables")
	  points(occAEA, pch = 19, cex = 0.5)

	  ## Panel 2: Cartoon of point versus spatial average climate estimate
	  par(mai = rep(0, 4))
	  
	    ### River 1
	    x1 <- seq(-1, 1, length.out = 1e2)
	    r1 <- function(x) exp(x) - exp(0) + sin(x * 10) / 10
	    y1 <- r1(x1)
	
	    ### River 2
	    x2 <- seq(0, 1, length.out = 1e2)
	    r2 <- function(x) exp(-x) - exp(0) + sin(x * 10) / 10
	    y2 <- r2(x2)
	
	  plot(0, 0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), axes = F)
	  clip(-1, 1, -1, 1)
	  points(x1, y1, type = "l", col = "blue", lwd = 5)
	  points(x2, y2, type = "l", col = "blue", lwd = 5)
	    
	  # Draw buffer area
	  x <- seq(-1, 0, length.out = 1e3)
	  clip(-2, 2, -2, 2)
	  polygon(c(x, rev(x)), c(sqrt(1 - x ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)
	  polygon(c(x + 1, rev(x + 1)), c(sqrt(1 - (x + 1) ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)
	  polygon(c(x, rev(x)), -c(sqrt(1 - x ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)
	  polygon(c(x + 1, rev(x + 1)), -c(sqrt(1 - (x + 1) ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)

	  # Draw point
	  points(0, 0, pch = 19, cex = 3)
	    
	  ## Panel 3: Climate-Latitude relationship
	  set.seed(01152016)
  	y <- (occ$Latitude - mean(occ$Latitude)) / sd(occ$Latitude) + 
		      rnorm(length(occ), 0, sqrt(0.1))
	  par(mai = c(0.75, 0.75, 0, 1/6))
	  plot(y, occ$Latitude, las = 1, axes = F, frame.plot = T, 
		  	 xlim = c(-2.5, 2.5), ylim = c(32, 45), cex.lab = 1.5,
			   pch = 19, col = "grey", xlab = "", ylab = "Latitude of origin")
	  title(xlab = "Climate variable", cex.lab = 1.5, line = 1)
	
  	for (i in seq(32, 44, 2)) {
		  axis(2, at = i, labels = bquote(.(i)*degree), las = 1)
	  	}
	  fit <- lm(occ$Latitude ~ y)
	  lines(range(y), range(predict(fit)), lwd = 2)

  # B. Example Climate-Trait variable, point estimated

  	## Panel 1: Map
	  par(mai = rep(0.5, 4))
	  plot(states, main = "B. Climate-Trait variables")
	  mtext(text = "point estimated", side = 3)
	  points(focAEA, pch = 19, cex = 0.5)

	  ## Panel 2: Cartoon of point versus spatial average climate estimate
	  par(mai = rep(0, 4))
	  plot(0, 0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), axes = F)
	  clip(-1, 1, -1, 1)
	  points(x1, y1, type = "l", col = "blue", lwd = 5)
	  points(x2, y2, type = "l", col = "blue", lwd = 5)
	  
	  # Draw buffer area
	  x <- seq(-1, 0, length.out = 1e3)
	  clip(-2, 2, -2, 2)
	  polygon(c(x, rev(x)), c(sqrt(1 - x ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)
	  polygon(c(x + 1, rev(x + 1)), c(sqrt(1 - (x + 1) ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)
	  polygon(c(x, rev(x)), -c(sqrt(1 - x ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)
	  polygon(c(x + 1, rev(x + 1)), -c(sqrt(1 - (x + 1) ^ 2), rep(1.1, 1e3)), 
	          col = "white", border = "white", lwd = 2)
	  
	  # Draw point
	  points(0, 0, pch = 19, cex = 3)
	  
	  ## Panel 3: Climate-trait relationship
	  par(mai = c(0.75, 0.75, 0, 1/6))
	  y <- (foc$Latitude - mean(foc$Latitude)) / sd(foc$Latitude) + 
		      rnorm(length(foc), 0, sqrt(0.1))
	  plot(y, foc$Latitude, las = 1, axes = F, frame.plot = T, 
		  	 xlim = c(-2.5, 2.5), ylim = c(32, 45),
			   pch = 19, col = "grey", xlab = "", ylab = "")
	  title(xlab = "Climate variable", ylab = "Trait", cex.lab = 1.5, line = 1)
	  fit <- lm(foc$Latitude ~ y)
	  lines(range(y), range(predict(fit)), lwd = 2)

  # C. Example Climate-Trait variable, spatially averaged

  	## Panel 1: Map
  	par(mai = rep(0.5, 4))
  	plot(states, main = "C. Climate-Trait variables")
  	tmp <- raster(states)
  	mtext(text = "spatially averaged", side = 3)
  	points(focAEA, pch = 21, cex = 2, col = "black", bg = "white")
	  points(focAEA, pch = 19, cex = 0.5)

  	## Panel 2: Cartoon of point versus spatial average climate estimate
	  par(mai = rep(0, 4))

  		### Climate grid

  			#### Functions to find closest point along "stream"
				fn <- function(par, x, y) {
					if (x < 0) {
						d1 <- sqrt((par - x) ^ 2 + (r1(par) - y) ^ 2)
						return(d1)
					} else {
						d1 <- sqrt((par - x) ^ 2 + (r1(par) - y) ^ 2)
						d2 <- sqrt((par - x) ^ 2 + (r2(par) - y) ^ 2)
						return(min(c(d1, d2)))
					}
				}

				getStreamDist <- function(x, y) {
					fit <- optim(x, fn, x = x, y = y, method = "Brent", lower = -2, 
											 upper = 2)
					retval <- setNames(c(fit$par, fit$value), c("s", "dist"))
					return(retval)
				}

			  #### "Climate" is fxn of distance from stream and pos along x-axis
			  xy <- expand.grid(x = seq(-0.95, 0.95, 0.1), y = seq(-0.95, 0.95, 0.1))
			  xyz <- cbind(xy, t(apply(xy, 1, function(X) getStreamDist(X[1], X[2]))))
			  xyz$z <- xyz$dist + xyz$x * 0.3
			  xyz$z <- xyz$z - min(xyz$z)
			  xyz$z <- xyz$z / max(xyz$z)

		  plot(0, 0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), axes = F)
		  clip(-1, 1, -1, 1)
		  points(x1, y1, type = "l", col = "blue", lwd = 5)
		  points(x2, y2, type = "l", col = "blue", lwd = 5)

		  apply(xyz, 1, function(X) rect(X["x"] - 0.05, X["y"] - 0.05, 
																	   X["x"] + 0.05, X["y"] + 0.05,
														 			   col = rgb(0, 0, 1, alpha = X["z"]), 
																	   border = NA))

			# Draw buffer area and circle
			x <- seq(-1, 0, length.out = 1e3)
			clip(-2, 2, -2, 2)
			polygon(c(x, rev(x)), c(sqrt(1 - x ^ 2), rep(1.1, 1e3)), 
			        col = "white", border = "white", lwd = 2)
			polygon(c(x + 1, rev(x + 1)), c(sqrt(1 - (x + 1) ^ 2), rep(1.1, 1e3)), 
			        col = "white", border = "white", lwd = 2)
			polygon(c(x, rev(x)), -c(sqrt(1 - x ^ 2), rep(1.1, 1e3)), 
			        col = "white", border = "white", lwd = 2)
			polygon(c(x + 1, rev(x + 1)), -c(sqrt(1 - (x + 1) ^ 2), rep(1.1, 1e3)), 
			        col = "white", border = "white", lwd = 2)
			
			x <- seq(-1, 1, length.out = 1e3)
			polygon(c(x, rev(x)), c(sqrt(1 - x ^ 2), -sqrt(1 - x ^ 2)),
			        border = "black", lwd = 2)
			
			# Draw point
      points(0, 0, pch = 19, cex = 3)
      
	## Panel 3: Climate-trait relationship

		par(mai = c(0.75, 0.75, 0, 1/6))
		plot(y, foc$Latitude, las = 1, axes = F, frame.plot = T, 
				 xlim = c(-2.5, 2.5), ylim = c(32, 45),
				 pch = 19, col = "grey", xlab = "", ylab = "")
		title(xlab = "Climate variable", ylab = "Trait", cex.lab = 1.5, line = 1)
		lines(range(y), range(predict(fit)), lwd = 2)

	dev.off()
