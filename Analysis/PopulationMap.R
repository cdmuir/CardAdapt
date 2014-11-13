#
# Niche modeling in Mimulus cardinalis
#

setwd("~/Dropbox/CardNiche")

# Records from GBIF (these already include UCJEPS, so maybe just use those?)

card <- gbif("Mimulus", "cardinalis")
card <- card[-which(card$lon > -80), ]
nrow(card)

# Jespon herbarium records
# Accessed July 3, 2013

card <- read.csv("~/Documents/Postdoc/CardNiche/card_jepson_records_minimal.csv", colClasses = c("factor", "integer", "integer", "factor", "character", "numeric", "numeric"))
with(card, plot(longitude, latitude))

# convert elevation strings to meters

x <- which(regexpr("[0-9]+ [-] [0-9]+ m", card$elevation) != -1)

le <- as.numeric(unlist(strsplit(card$elevation[x], " - "))[seq(1, 2 * length(x), 2)])
ue <- as.numeric(unlist(strsplit(unlist(strsplit(card$elevation[x], " - "))[seq(2, 2 * length(x), 2)], " m")))
card$elevation[x] <- (le + ue) / 2

x <- which(regexpr("[0-9]+ [-] [0-9]+ ft", card$elevation) != -1)
le <- as.numeric(unlist(strsplit(card$elevation[x], " - "))[seq(1, 2 * length(x), 2)])
ue <- as.numeric(unlist(strsplit(unlist(strsplit(card$elevation[x], " - "))[seq(2, 2 * length(x), 2)], " ft")))
card$elevation[x] <- (le + ue) / 2

x <- regexpr(" ft", card$elevation)

card$elevation[x != -1] <- as.numeric(substr(card$elevation[x != -1], 1, x[x != -1] - 1)) * 0.3048

x <- regexpr(" m", card$elevation)

card$elevation[x != -1] <- as.numeric(substr(card$elevation[x != -1], 1, x[x != -1] - 1))

card$elevation <- as.numeric(card$elevation)

#
## Anget lab seed collections
#

library(raster)
library(dismo)

card2 <- read.csv("~/Dropbox/CardNiche/CardPopulations.csv")

# get rid of distant pop
card2 <- card2[card2$Lon < -115 & !is.na(card2$Lon), ]
sp2 <- SpatialPoints(cbind(card2$Lon, card2$Lat))

#
## BIOCLIM
#

# files <- list.files("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/bio_30s_bil", pattern = "bil")
# bioclim <- stack(paste("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/bio_30s_bil/", files, sep = ""))

# files <- list.files("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/tmin_30s_bil", pattern = "bil")
# tmin <- stack(paste("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/tmin_30s_bil/", files, sep = ""))

# files <- list.files("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/tmax_30s_bil", pattern = "bil")
# tmax <- stack(paste("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/tmax_30s_bil/", files, sep = ""))

# files <- list.files("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/tmean_30s_bil", pattern = "bil")
# tmean <- stack(paste("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/tmean_30s_bil/", files, sep = ""))

# files <- list.files("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/prec_30s_bil", pattern = "bil")
# prec <- stack(paste("~/Documents/GradSchool/Projects/Solanum_ENM/PredictorData/prec_30s_bil/", files, sep = ""))

#
# Extract bioclim
#

# bioclim2 <- data.frame(extract(bioclim, sp2))
# tmin.by.month <- data.frame(extract(tmin, sp2))
# tmax.by.month <- data.frame(extract(tmax, sp2))
# tmean.by.month <- data.frame(extract(tmean, sp2))
# prec.by.month <- data.frame(extract(prec, sp2))

# card2 <- cbind(card2, bioclim2, tmin.by.month, tmax.by.month, tmean.by.month, prec.by.month)

#
## Write new csv with bioclim data
#

# write.csv(card2, file = "CardPopulations_with_bioclim.csv", row.names = F)
card2 <- read.csv("CardPopulations_with_bioclim.csv")

#
## Plot focal populations
#

load("USA_adm1.RData")

calore <- gadm[gadm$NAME_1 == "California" | gadm$NAME_1 == "Oregon", ]

pdf("FocalPopulations.pdf", 5, 5)
plot(calore)

focal <- read.csv("PopulationSelection.csv")
levels(focal$Group)
focal$Group <- factor(as.character(focal$Group), levels = c("South Margin", "Transverse", "South Sierras", "Central Sierras", "North Sierras", "North Coast", "North Margin"))

palette(heat.colors(7))
with(focal, points(Lon, Lat, pch = 21, col = "black", bg = Group))
dev.off()

with(card2, plot(tmax_7, prec_7))
with(focal, points(tmax_7, prec_7, pch = 21, col = "black", bg = Group))

# fit seasonal temperature function (gaussian)

x <- numeric(nrow(tmean.by.month))
out <- data.frame(a = x, u = x, s = x, h = x, r2 = x)
for (i in 1:nrow(tmean.by.month))
{
	data <- data.frame(x = c(1, 10:12, 2:9), y = unlist(tmean.by.month[i, ]))
	fit <- nls(y ~ a * exp(-(x - u) ^ 2 / (2 * s ^ 2)) + h, 
		start = list(a = diff(range(data$y)), u = 6, s = 2, h = min(data$y)), data = data)
	out[i, ] <- c(coef(fit), cor(data$y, predict(fit)) ^ 2)
}

# Proportion of year above 20 degrees

x <- seq(1, 12, length.out = 100)
out$p20 <- numeric(nrow(out))
for (i in 1:nrow(out))
{
	out$p20[i] <- with(out[i, ], length(which(a * exp(-(x - u) ^ 2 / (2 * s ^ 2)) + h > 200)) / length(x))
}

with(card2, plot(bio_7, bio_17))

with(card2[card2$Demography.site. == "Yes", ], points(bio_7, bio_17, pch = 19, col = "red"))
with(card2[card2$Population.genomic.dataset. == "Yes", ], points(bio_7, bio_17, pch = 19, col = "blue"))
with(card2[card2$Demography.site. == "Yes" & card2$Population.genomic.dataset. == "Yes", ], points(bio_7, bio_17, pch = 19, col = "purple"))

plot(card2$tmean_7, card2$prec_7, log = "")
with(card2[card2$Demography.site. == "Yes", ], points(tmean_7, prec_7, pch = 19, col = "red"))
with(card2[card2$Demography.site. == "Yes", ], text(tmean_7, prec_7, labels = X))
with(card2[card2$Population.genomic.dataset. == "Yes", ], points(tmean_7, prec_7, pch = 19, col = "blue"))
with(card2[card2$Population.genomic.dataset. == "Yes", ], text(tmean_7, prec_7, labels = X))

with(card2[card2$Demography.site. == "Yes" & card2$Population.genomic.dataset. == "Yes", ], points(tmean_7, prec_7, pch = 19, col = "purple"))

card2$X[card2$Demography.site. == "Yes" & card2$Population.genomic.dataset. == "Yes"]

with(card2, plot(Lon, Lat))
plot(log(bioclim[[10]] + 1), add = T)

with(card2[card2$Demography.site. == "Yes", ], points(Lon, Lat, pch = 19, col = "red"))
with(card2[card2$Population.genomic.dataset. == "Yes", ], points(Lon, Lat, pch = 19, col = "blue"))
with(card2[card2$Demography.site. == "Yes" & card2$Population.genomic.dataset. == "Yes", ], points(Lon, Lat, pch = 19, col = "purple"))


#
## Find crosses of focal accessions
#

card1 <- read.csv("~/Documents/Postdoc/CardNiche/Angert_Lab_Mimulus_Accessions_2012_working_version.csv", nrows = 2334)

crosses <- read.csv("~/Documents/Postdoc/CardNiche/Angert_Lab_Mimulus_Accessions_2012_working_version.csv", skip = 2334, nrows = 2165, header = F)

focal <- read.csv("~/Documents/Postdoc/CardNiche/FocalPopulations.csv")

sapply(focal$Accession, function(X) length(grep(X, crosses[, 12])))



plot(card2$Lat, bioclim2$bio_9) # mean temp of driest quarter
plot(card2$Lat[bioclim2$alt < 1000], tmean.by.month$tmean_8[bioclim2$alt < 1000])
plot(card2$Lat[bioclim2$alt < 1000], tmax.by.month$tmax_7[bioclim2$alt < 1000])

plot(card2$Lat[bioclim2$alt < 1000], bioclim2$bio_9[bioclim2$alt < 1000])
plot(card2$Lat, bioclim2$bio_1)

# Prin comp of bioclim
prbio <- princomp(bioclim2, cor = T)
summary(prbio)
plot(bioclim2$alt, card2$Elev)

plot(1:12, prec.by.month[100, c(1, 5:12, 2:4)], type = "l", ylim = c(0, 250))
points(1:12, tmin.by.month[2, c(1, 5:12, 2:4)], type = "l", lty = 2)

# prec by month
prec.by.month <- data.frame(extract(prec, sp2[c(1, 2330), ]))
plot(1:12, prec.by.month[1, c(1, 5:12, 2:4)], type = "l", ylim = c(0, 150))
points(1:12, prec.by.month[2, c(1, 5:12, 2:4)], type = "l", lty = 2)

with(card2, table(Demography.site., Population.genomic.dataset.))



plot(card$elevation, bioclim[, 1])
plot(card2$Elevation.Converted, bioclim2[, 1])

# PLOT

plot(predictors[[15]], xlim = range(card$longitude, na.rm = T), ylim = range(card2$Latitide, na.rm = T))

with(card, points(longitude, latitude, pch = 19))
with(card2, points(Longitude, Latitide, pch = 19, col = "red"))

sort(tapply(card2$Latitide, card2$Site, mean, na.rm = T))