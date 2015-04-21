setwd("~/Google Drive/CardLocalAdaptation/Data/Climate/ClimateWNA")
library(raster)
# environmental
# growing degree days, summer precip, drought indicies
# variation within and between sites

#make big stack
files <- list.files("NORM_6190_Bioclim_ASCII", pattern = ".asc")
bioclim <- stack(paste("NORM_6190_Bioclim_ASCII/", files, sep = "")[1])

latlon <- SpatialPoints(popenv[, c("Lon", "Lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))


popenv$ggd <- 
tmp <- raster::extract(bioclim, latlon[1, ], buffer = 1e6, cellnumbers = T) # buffer is in meters

r <- bioclim[[1]]
r[] <- NA
r[tmp[[1]][, 1]] <- 1
plot(bioclim)

with(popenv, plot(Lat, ggd))
	plot(popenv$ggd, betas$lsmeans.table$Estimate[x], log = "")
	
cor(cbind(popenv[, 75], betas$lsmeans.table$Estimate[x]))

# get all card points to get average selective environment along with latitude.


# PRISM file name format:
PRISM_ppt_stable_4kmD1_YYYYMMDD_bil
PRISM_tmax_stable_4kmD1_YYYYMMDD_bil
PRISM_tmin_stable_4kmD1_YYYYMMDD_bil
PRISM_tmean_stable_4kmD1_YYYYMMDD_bil

setwd("~/Google Drive/CardLocalAdaptation/Data/Climate/PRISM_daily")
files <- list.files(pattern = ".zip")[1:2]
for (i in files)
{
	system(paste("unzip", i))
}