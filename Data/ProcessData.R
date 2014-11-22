setwd("~/Google Drive/CardLocalAdaptation/Data")

library(lme4)
library(reshape2)
library(lubridate)
library(MASS)
	
#
#	Preparing datasheet for LLL analysis
#

	data <- read.csv("MasterDatasheet_in.csv", nrows = 768)
	cols2keep <- c("indiv", "Block", "Col", "Row", "Population", "Family", "TempTrt",
		"WaterTrt", "LLL_0512", "LLL_0515", "LLL_0520", "LLL_0523", "LLL_0526", 
		"LLL_0529", "LLL_0602", "LLL_0605", "LLL_0609", "LLL_0612")
	data$Row <- as.character(data$Row)
	data$Population[data$Population == "Sweetwater River"] <- "Cuyamaca Rancho"
	data$Population <- factor(as.character(data$Population))
	data$indiv <- gl(768, 1)
	data <- data[data$UseForLLLGrowth, cols2keep]
	data$TempTrt <- factor(as.character(data$TempTrt))
	data$WaterTrt <- factor(as.character(data$WaterTrt))

	# Remove all < 1 mm, call NA
	data$LLL_0512 <- as.numeric(gsub("< 1", "NA", data$LLL_0512))
	data$LLL_0515 <- as.numeric(gsub("< 1", "NA", data$LLL_0515))
	data$LLL_0520 <- as.numeric(gsub("< 1", "NA", data$LLL_0520))
	data$LLL_0523 <- as.numeric(gsub("< 1", "NA", data$LLL_0523))
	data$LLL_0526 <- as.numeric(gsub("< 1", "NA", data$LLL_0526))

#
#	LLL growth rates
#

	# Melt data.frame
	lll.data <- melt(data, value.name = "LLL")
	lll.data$Day <- mdy(paste(substr(lll.data$variable, 5, 8), "2014", sep = ""))
	lll.data$DayN <- (as.numeric(lll.data$Day) - 1399852800) / 86400
	
	# Get rid of NAs
	lll.data <- lll.data[!is.na(lll.data$LLL), ]
		
	# Refactor indiv
	lll.data$indiv <- as.factor(as.character(lll.data$indiv))

	#
	# Empirical Bayes Approach: get coefficients for every plant
	#

	mm <- lmer(log(LLL) ~ poly(DayN, 2, raw = T) + (poly(DayN, 2, raw = T)|indiv),
		data = lll.data)
	# Logistic not a very good fit. Something else worth trying?
	# mm <- nlmer(log(LLL) ~ SSlogis(DayN, Asym, xmid, scal) ~ Asym|indiv,
		# data = lll.data, start = c(Asym = 5, xmid = 15, scal = 10))
	
	# TEMP: check for outliers
	plot(mm) # no outliers, but residuals are odd.

	# Extract coefficients and add to master data list
	Y <- t(unlist(fixef(mm)) + t(ranef(mm)$indiv))[match(data$indiv, 
		rownames(ranef(mm)$indiv)), ]
	colnames(Y) <- c("b0_lll", "b1_lll", "b2_lll")
	data <- cbind(data, Y)
		
	# For visual representation, plot LLL vs DayN with curves
	pdf("LLL by individual.pdf", 4, 4)
	r2 <- numeric(nlevels(lll.data$indiv))
	names(r2) <- levels(lll.data$indiv)
	par(mar = c(5, 5, 1, 1))
	for (i in levels(lll.data$indiv))
	{
		plot(0, 0, xlim = c(0, 31), ylim = c(0, 5), type = "n", las = 1,
			xlab = "Day", ylab = "log(LLL)", cex.lab = 1.5)
		with(subset(lll.data, lll.data$indiv == i), points(DayN, log(LLL), pch = 19))
		with(subset(lll.data, lll.data$indiv == i), points(DayN, 
			predict(mm)[lll.data$indiv == i], type = "l", col = "blue"))
		mtext(text = paste("indiv:", i))
		r2[i] <- with(subset(lll.data, lll.data$indiv == i), cor(log(LLL), 
			predict(mm)[lll.data$indiv == i]))
		mtext(text = bquote(r^2 == .(r2[i])), line = -1)
	}
	dev.off()

#
#	Preparing datasheet for Height analysis
#

	data2 <- read.csv("MasterDatasheet_in.csv", nrows = 768)
	cols2keep <- c("indiv", "Population", "Family", "Transplant", "TempTrt", "WaterTrt",
		"Height_0529", "Height_0602", "Height_0605", "Height_0609", "Height_0612", 
		"Height_0616", "Height_0620")
	data2$Population[data2$Population == "Sweetwater River"] <- "Cuyamaca Rancho"
	data2$Population <- factor(as.character(data2$Population))
	data2$indiv <- gl(768, 1)
	data2 <- data2[data2$UseForHeightGrowth, cols2keep]
	data2$TempTrt <- factor(as.character(data2$TempTrt))
	data2$WaterTrt <- factor(as.character(data2$WaterTrt))

	# Melt data.frame
	height.data <- melt(data2)
	colnames(height.data)[colnames(height.data) == "value"]  <- "height"
	height.data$Day <- mdy(paste(substr(height.data$variable, 8, 11), "2014", sep = ""))
	height.data$DayN <- (as.numeric(height.data$Day) - 1399852800) / 86400
	
	# Get rid of NAs
	height.data <- height.data[!is.na(height.data$height), ]
	
	# Remove those with less than 3 datapoints (get rid of this)
	height.data <- height.data[!(height.data$indiv %in% 
		which(table(height.data$indiv) < 3)), ]
	
	# Refactor indiv
	height.data$indiv <- as.factor(as.character(height.data$indiv))

	
	#
	# Empirical Bayes Approach: get coefficients for every plant
	#

	mm <- lmer(log(height + 1) ~ poly(DayN, 2, raw = T) + (poly(DayN, 2, raw = T)|indiv),
		data = height.data)
	
	# Extract coefficients and add to master data list
	Y <- t(unlist(fixef(mm)) + t(ranef(mm)$indiv))[match(data2$indiv, 
		rownames(ranef(mm)$indiv)), ]
	colnames(Y) <- c("b0_height", "b1_height", "b2_height")
	data2 <- cbind(data2, Y)

	pdf("height by individual.pdf", 4, 4)
	r2 <- numeric(nlevels(height.data$indiv))
	names(r2) <- levels(height.data $indiv)
	par(mar = c(5, 5, 1, 1))
	for (i in levels(height.data$indiv))
	{
		plot(0, 0, xlim = c(17, 39), ylim = c(0, 4), type = "n", las = 1,
			xlab = "Day", ylab = "height", cex.lab = 1.5)
		with(subset(height.data, height.data$indiv == i), points(DayN, log(height + 1),
			pch = 19))
		with(subset(height.data, height.data$indiv == i), points(DayN, 
			predict(mm)[height.data$indiv == i], type = "l", col = "blue"))
		mtext(text = paste("indiv:", i))
		r2[i] <- with(subset(height.data, height.data $indiv == i), cor(log(height + 1), 
			predict(mm)[height.data$indiv == i]))
		mtext(text = bquote(r^2 == .(r2[i])), line = -1)
	}
	dev.off()

	#
	#	Combine lll and height data, output as MasterDatasheet_out.csv
	#

	data$Height_0529 <- data2$Height_0529[match(data$indiv, data2$indiv)]
	data$Height_0602 <- data2$Height_0602[match(data$indiv, data2$indiv)]
	data$Height_0605 <- data2$Height_0605[match(data$indiv, data2$indiv)]
	data$Height_0609 <- data2$Height_0609[match(data$indiv, data2$indiv)]
	data$Height_0612 <- data2$Height_0612[match(data$indiv, data2$indiv)]
	data$Height_0616 <- data2$Height_0616[match(data$indiv, data2$indiv)]
	data$Height_0620 <- data2$Height_0620[match(data$indiv, data2$indiv)]
	data$b0_height <- data2$b0_height[match(data$indiv, data2$indiv)]
	data$b1_height <- data2$b1_height[match(data$indiv, data2$indiv)]
	data$b2_height <- data2$b2_height[match(data$indiv, data2$indiv)]

	# Model-based lll growth
	# Cool: Day 0 - 31
	# Hot: Day 0 - 24
	data$lll_start <- exp(data$b0_lll)
	data$lll_end <- exp(data$b0_lll + 
					data$b1_lll * ifelse(data$TempTrt == "Cool", 31, 24) + 
					data$b2_lll * ifelse(data$TempTrt == "Cool", 31 ^ 2, 24 ^ 2))
	data$lll_AbsGrowth <- (data$lll_end - data$lll_start) / 
		ifelse(data$TempTrt == "Cool", 31, 24)
	data$lll_RelGrowth <- (log(data$lll_end) - log(data$lll_start)) / 
		ifelse(data$TempTrt == "Cool", 31, 24)

	# Model-based height growth
	# Cool: Day 17 - 39
	# Hot: Day 17 - 31
	data$height_start <- exp(data$b0_height + data$b1_height * 17 + 
		data$b2_height * 17 ^ 2) - 1
	data$height_end <- exp(data$b0_height + 
		data$b1_height * ifelse(data$TempTrt == "Cool", 39, 31) + 
		data$b2_height * ifelse(data$TempTrt == "Cool", 39 ^ 2, 31 ^ 2)) - 1
	data$height_AbsGrowth <- (data$height_end - data$height_start) / 
		ifelse(data$TempTrt == "Cool", 22, 14)
	data$height_RelGrowth <- (log(data$height_end) - log(data$height_start)) / 
		ifelse(data$TempTrt == "Cool", 22, 14)

#
#	Add photosynthetic data
#

	licor <- read.csv("CardLICOR.csv", skip = 1)
	data$ID <- with(data, paste(Block, Col, Row))
	data$ID <- as.character(data$ID)
	licor$Block <- with(licor, paste(ifelse(Tray <= 4, "MM5-", "MM7-"), Tray, sep = ""))
	licor$ID <- with(licor, paste(Block, Col, Row))
	licor$ID <- as.character(licor$ID)

	# Average Photo and Cond by ID
	avgPhoto <- tapply(licor$Photo, licor$ID, mean)
	avgCond <- tapply(licor$Cond, licor$ID, mean)
	avgTleaf <- tapply(licor$Tleaf, licor$ID, mean)
	
	# Remove photo data for plants not in master data
	# (these were transplants, which have been excluded)
	avgPhoto <- avgPhoto[names(avgPhoto) %in% data$ID]
	avgCond <- avgCond[names(avgCond) %in% data$ID]
	avgTleaf <- avgTleaf[names(avgTleaf) %in% data$ID]

	data$Photo <- data$Cond <- data$Tleaf <- rep(NA, nrow(data))
	data$Photo[na.omit(match(names(avgPhoto), data$ID))] <- avgPhoto
	data$Cond[match(names(avgCond), data$ID)] <- avgCond
	data$Tleaf[match(names(avgTleaf), data$ID)] <- avgTleaf

	# 'Residual' photosynthetic rate after accounting for conductance
	fit <- stepAIC(lm(log(Photo) ~ log(Cond) * TempTrt * WaterTrt,
		data = data[data$Photo > 0 & !is.na(data$Photo), ]))
	data$resPhoto <- log(data$Photo) - predict(fit, 
		new = data[, c("Cond", "TempTrt", "WaterTrt")])	

#
#	Add additional columns from MasterDatasheet_in
#

	tmp <- read.csv("MasterDatasheet_in.csv", nrows = 768)
	tmp$Row <- as.character(tmp$Row)
	tmp$Population[tmp$Population == "Sweetwater River"] <- "Cuyamaca Rancho"
	tmp$Population <- factor(as.character(tmp$Population))
	tmp$indiv <- gl(768, 1)
	tmp$TempTrt <- factor(as.character(tmp$TempTrt))
	tmp$WaterTrt <- factor(as.character(tmp$WaterTrt))

	tmp <- subset(tmp, tmp$UseForLLLGrowth)
	# Should be true if tmp and data are in same order
	all(tmp$indiv == data$indiv)

	# add biomass and flowers at harvest as well as harvest date [STILL NEEDS QC]
	data$LeavesDW_g <- ifelse(tmp$UseForBiomass, tmp$LeavesDW_g, NA)
	data$ShootsDW_g <- ifelse(tmp$UseForBiomass, tmp$ShootsDW_g, NA)
	data$RootsDW_g <- ifelse(tmp$UseForBiomass, tmp$RootsDW_g, NA)
	data$Biomass <- with(data, LeavesDW_g + ShootsDW_g + RootsDW_g)
	data$NFlwr <- ifelse(tmp$UseForBiomass, tmp$NFlwr, NA)
	data$NBud <- ifelse(tmp$UseForBiomass, tmp$NBud, NA)
	data$HarvestDate <- ifelse(nchar(as.character(tmp$HarvestDate)) == 6, 
		(as.numeric(dmy(paste(as.character(tmp$HarvestDate), "-2014", sep = ""))) - 
		1399852800) / 86400, NA)
	
	# add germination date
	data$MinGermDay <- tmp$MinGermDay
	data$MaxGermDay <- tmp$MaxGermDay
	data$AvgGermDay <- (tmp$MinGermDay + tmp$MaxGermDay) / 2

	# add Mortality
	data$DiedFromStress <- tmp$DiedFromStress
	
#
#	Export for Master Analysis
#

	data <- data[which(!(is.na(data$lll_AbsGrowth) | is.na(data$height_AbsGrowth))), ]
	write.csv(data, "MasterDatasheet_out.csv")

	
### stuff to move possibly	
	#
	#	start fitting with LMER
	#

	# clean up/write function/make sure proper tests are being done (double check Bolker et al. 2009)
	fit.b0.1 <- lmer(b0_eb ~ TempTrt * WaterTrt * Population + (1|Family), 
		data = subset(data, !is.na(data$b0_eb)))
	# population interaction effects
	fit.b0.2 <- update(fit.b0.1, . ~ . - TempTrt:WaterTrt:Population)
	KRmodcomp(fit.b0.1, fit.b0.2)
	fit.b0.3 <- update(fit.b0.2, . ~ . - TempTrt:Population)
	KRmodcomp(fit.b0.2, fit.b0.3)
	fit.b0.4 <- update(fit.b0.2, . ~ . - WaterTrt:Population)
	KRmodcomp(fit.b0.2, fit.b0.4)

	# Population effect
	fit.b0.5 <- lmer(b0_eb ~ TempTrt * WaterTrt + Population + (1|Family), 
		data = subset(data, !is.na(data$b0_eb))) # this is best model, I think
	fit.b0.6 <- update(fit.b0.5, . ~ . - Population)
	KRmodcomp(fit.b0.5, fit.b0.6) # keep population as factor

	# Treatment effects
	fit.b0.7 <- update(fit.b0.5, . ~ . - TempTrt:WaterTrt)
	KRmodcomp(fit.b0.5, fit.b0.7) # keep temp x water interaction

	fit.b0.8 <- update(fit.b0.5, . ~ . - WaterTrt - TempTrt:WaterTrt)
	KRmodcomp(fit.b0.5, fit.b0.8) # keep water

	fit.b0.9 <- update(fit.b0.5, . ~ . - TempTrt - TempTrt:WaterTrt)
	KRmodcomp(fit.b0.5, fit.b0.9) # keep temp interaction

	# now do with b1
	fit.b1.1 <- lmer(b1_eb ~ TempTrt * WaterTrt * Population + (1|Family), 
		data = subset(data, !is.na(data$b0_eb)))
	# population interaction effects
	fit.b1.2 <- update(fit.b1.1, . ~ . - TempTrt:WaterTrt:Population)
	KRmodcomp(fit.b1.1, fit.b1.2)
	fit.b1.3 <- update(fit.b1.2, . ~ . - TempTrt:Population)
	KRmodcomp(fit.b1.2, fit.b1.3) # marginal sig
	fit.b1.4 <- update(fit.b1.2, . ~ . - WaterTrt:Population)
	KRmodcomp(fit.b1.2, fit.b1.4)

	# Population effect
	fit.b1.5 <- lmer(b1_eb ~ TempTrt * WaterTrt + Population + (1|Family), 
		data = subset(data, !is.na(data$b1_eb))) # this is best model, I think
	fit.b1.6 <- update(fit.b1.5, . ~ . - Population)
	KRmodcomp(fit.b1.5, fit.b1.6) # keep population as factor

	# Treatment effects
	fit.b1.7 <- update(fit.b1.5, . ~ . - TempTrt:WaterTrt)
	KRmodcomp(fit.b1.5, fit.b1.7) # keep temp x water interaction

	fit.b1.8 <- update(fit.b1.5, . ~ . - WaterTrt - TempTrt:WaterTrt)
	KRmodcomp(fit.b1.5, fit.b1.8) # keep water

	fit.b1.9 <- update(fit.b1.5, . ~ . - TempTrt - TempTrt:WaterTrt)
	KRmodcomp(fit.b1.5, fit.b1.9) # keep temp

	# now do with b2
	fit.b2.1 <- lmer(b2_eb ~ TempTrt * WaterTrt * Population + (1|Family), 
		data = subset(data, !is.na(data$b0_eb)))
	# population interaction effects
	fit.b2.2 <- update(fit.b2.1, . ~ . - TempTrt:WaterTrt:Population)
	KRmodcomp(fit.b2.1, fit.b2.2)
	fit.b2.3 <- update(fit.b2.2, . ~ . - TempTrt:Population)
	KRmodcomp(fit.b2.2, fit.b2.3) # highly significant
	fit.b2.4 <- update(fit.b2.2, . ~ . - WaterTrt:Population)
	KRmodcomp(fit.b2.2, fit.b2.4)

	# Population effect
	fit.b2.5 <- lmer(b2_eb ~ TempTrt * WaterTrt + Population + (1|Family), 
		data = subset(data, !is.na(data$b2_eb))) # this is best model, I think
	fit.b2.6 <- update(fit.b2.5, . ~ . - Population)
	KRmodcomp(fit.b2.5, fit.b2.6) # keep population as factor

	# Treatment effects
	fit.b2.7 <- update(fit.b2.5, . ~ . - TempTrt:WaterTrt)
	KRmodcomp(fit.b2.5, fit.b2.7) # keep temp x water interaction

	fit.b2.8 <- update(fit.b2.5, . ~ . - WaterTrt - TempTrt:WaterTrt)
	KRmodcomp(fit.b2.5, fit.b2.8) # keep water

	fit.b2.9 <- update(fit.b2.5, . ~ . - TempTrt - TempTrt:WaterTrt)
	KRmodcomp(fit.b2.5, fit.b2.9) # keep temp

	#
	#	Latitude versus b0, b1, b2
	#
	
	fit.b0 <- lmer(b0_height ~ -1 + Population + TempTrt * WaterTrt + (1|Family), 
		data = subset(data, !is.na(data$b2_height)))
	fit.b1 <- lmer(b1_height ~ -1 + Population + TempTrt * WaterTrt + (1|Family), 
		data = subset(data, !is.na(data$b1_height)))
	fit.b2 <- lmer(b2_height ~ -1 + Population + TempTrt * WaterTrt + (1|Family), 
		data = subset(data, !is.na(data$b2_height)))
	
	popenv$height_b0 <- fixef(fit.b0)[1:16][match(popenv$Site, 
		substr(names(fixef(fit.b0)[1:16]), 11, nchar(names(fixef(fit.b0)[1:16]))))]
	popenv$height_b1 <- fixef(fit.b1)[1:16][match(popenv$Site, 
		substr(names(fixef(fit.b1)[1:16]), 11, nchar(names(fixef(fit.b1)[1:16]))))]
	popenv$height_b2 <- fixef(fit.b2)[1:16][match(popenv$Site, 
		substr(names(fixef(fit.b2)[1:16]), 11, nchar(names(fixef(fit.b2)[1:16]))))]
	
	
	# MCMCglmm (TBD)
	prior1 <- list(R = list(V = diag(3), n = 3), G = list(G1 = list(V = diag(2), n = 3), G2 = list(V = diag(2), n = 3))) 

	mm1 <- MCMCglmm(cbind(b0_eb, b1_eb, b2_eb) ~ TempTrt * WaterTrt, random = ~ us(1 + TempTrt):Family + us(1 + TempTrt):Population, rcov = ~ us(trait):units, prior = prior1, data = data, family = rep("gaussian", 3))


#
#	put height and LLL by individual on same axis
#

	x <- rnorm(1e3, -2, 2); y <- rnorm(1e3, 2, 4)
	rescale <- function(x, y)
	{
		ma <- mean(y) - mean(x)
		va <- sd(y) / sd(x)
		y1 <- (y - ma) / va
		return(y1)
	}

	pdf("lll+height by individual.pdf", 4, 4)
	par(mar = c(5, 5, 1, 1))
	
	height.data$reheight <- rescale(log(lll.data$LLL), log(height.data$height + 1))
	for (i in intersect(lll.data$indiv, height.data$indiv))
	{
		with(subset(lll.data, lll.data$indiv == i), plot(DayN, log(LLL), pch = 19,
			col = "blue", xlim = c(0, 39), ylim = c(0, 5), las = 1, main = i))
		with(subset(height.data, height.data$indiv == i), points(DayN, log(height + 1),
			pch = 19, col = "red", xlim = c(0, 39), ylim = c(0, 5), las = 1))
	}
	dev.off()


plot(fit)

palette(grey(seq(0, 1, length.out = 16)))
pdf("LLL_by_pop.pdf", 10, 5)
par(mar = c(5, 5, 1, 1), mfrow = c(1, 2))

for (i in 2:18)
{
	pop <- levels(data$Population)[i]
	with(subset(data, data$Population == pop & data$TempTrt == "Hot"),
		plot(DayN, LLL, type = "p", pch = 19, log = "", main = pop, 
		xlim = c(min(data$DayN), max(data$DayN)), ylim = c(min(data$LLL, na.rm = T), 
		max(data$LLL, na.rm = T))))

	with(subset(data, data$Population == pop & data$TempTrt == "Cool"),
		plot(DayN, LLL, type = "p", pch = 19, log = "", main = pop, 
		xlim = c(min(data$DayN), max(data$DayN)), ylim = c(min(data$LLL, na.rm = T), 
		max(data$LLL, na.rm = T))))
}

dev.off()	
	
# Height versus LLL:

with(data, plot(LLL_0529, Height_0529 + 1, log = "", pch = 19, col = TempTrt))

with(data[600, ], plot(c(12, 15, 20, 23, 26, 29), c(LLL_1, LLL_2, LLL_3, LLL_4, LLL_5, LLL_6), type = "b", log = "y"))

summary(fit <- aov(log(Growth) ~ Block * Population/Family, data = data[!is.na(data$Growth), ]))
summary(fit <- boxplot(log(LLL_2) ~ Block, data = data[!is.na(data$LLL_2), ]))
boxplot(log(Growth) ~ Block, data = data[!is.na(data$LLL_2), ])
y <- data$MinGermDay[!is.na(data$MinGermDay)]
fit <- glm(y ~ Block + Population, data = data[!is.na(data$MinGermDay), ], family = "poisson")

summary(mm1 <- lmer(log(Growth) ~ 1 + (1|Population), data = data[!is.na(data$LLL_1), ]))

mm2 <- MCMCglmm(log(Growth) ~ Block, random = ~ Population, family = "gaussian", data = data)


mm2 <- MCMCglmm(cbind(MinGermDay, MaxGermDay) ~ 1, random = ~ Block + Family + Population : Family, family = "cenpoisson", data = data)






popenv <- popenv[order(popenv$Site), ]

lll <- tapply(data$LLL_0523, data$Population, mean, na.rm = T)
lll <- tapply(data$LLL_0523, data$Population, sd, na.rm = T)
	plot(popenv$Lat, lll[match(popenv$Site, rownames(x0))])
	cor.test(popenv$Lat, lll[match(popenv$Site, rownames(x0))])

fit <- lm(log(LLL_0523) ~ -1 + Population + Population:TempTrt + Population:WaterTrt,
	data = data)

# FIGURES FOR KILLAM APPLICATION
pdf("GrowthChamberData.pdf", 4, 4)
par(mar = c(5, 5, 1, 1))
plot(popenv$Lat, coef(fit)[1:16], xlim = c(32, 44), ylim = c(1.8, 3.1), pch = 19,
	las = 1, cex.lab = 1.5, xlab = "Latitude", ylab = "log(Longest leaf length)")
for (i in 1:16) arrows(popenv$Lat[i], confint(fit)[i, 1], popenv$Lat[i], 
	confint(fit)[i, 2], length = 0.05, angle = 90, code = 3)
abline(lm(coef(fit)[1:16] ~ popenv$Lat))

plot(popenv$Lat, coef(fit)[17:32], xlim = c(32, 44), ylim = c(-0.1, 1), pch = 19,
	las = 1, cex.lab = 1.5, xlab = "Latitude", ylab = "Temperature Response")
for (i in 1:16) arrows(popenv$Lat[i], confint(fit)[i + 16, 1], popenv$Lat[i], 
	confint(fit)[i + 16, 2], length = 0.05, angle = 90, code = 3)
abline(lm(coef(fit)[17:32] ~ 1))

plot(popenv$Lat, coef(fit)[33:48], xlim = c(32, 44), ylim = c(-1, 0.5), pch = 19,
	las = 1, cex.lab = 1.5, xlab = "Latitude", ylab = "Soil Moisture Response")
for (i in 1:16) arrows(popenv$Lat[i], confint(fit)[i + 32, 1], popenv$Lat[i], 
	confint(fit)[i + 32, 2], length = 0.05, angle = 90, code = 3)
abline(lm(coef(fit)[33:48] ~ 1))

dev.off()

fit <- lm(log(Height_0609 + 1) ~ -1 + Population + Population:TempTrt + Population:WaterTrt, data = data2)

plot(popenv$Lat, coef(fit)[1:16])
plot(popenv$Lat, coef(fit)[17:32])
plot(popenv$Lat, coef(fit)[33:48])

height <- tapply(log(data2$Height_0609 + 1), data2$Population, mean, na.rm = T)
	plot(popenv$Lat, height[match(popenv$Site, rownames(x0))])
	cor.test(popenv$Lat, height[match(popenv$Site, rownames(x0))])
