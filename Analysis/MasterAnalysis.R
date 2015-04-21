#
#	Master Analysis of M. cardinalis Growth Chamber Experiment (Spring 2014)
#
#	Title:	 Local Physiological adaptaiton in Mimulus cardinalis
#	Authors: Christopher D. Muir and Amy L. Angert
#	Journal: Evolution?
#


## PART 1:	TRAIT VARIATION
## PART 2: 	SELECTIVE AGENT
## PART 3	???

### Preliminaries

#
##	Libraries
#

	library(lme4)
	library(MCMCglmm)
	library(MASS)
	library(car)
	library(pbkrtest)
	library(lmerTest)
	library(survival)
	library(VSURF)
	
#
##	Set working directory
#

	setwd("~/Google Drive/CardLocalAdaptation")

#
##	Source custom functions
#

	source("Analysis/functions.R")
	
#
##	Focal population climate data
#
	# This data.frame will change as climate data pipeline is updated
	# Populations are ordered alphabetically by site to match output from stats
	load('~/Google Drive/CardLocalAdaptation/Data/Climate/Robjects/focClim.RData')
	
#
##	Read in processed datasets
#

	data <- read.csv("Data/MasterDatasheet_out.csv")

#
##	Difference in germiantion timing (yes) and intitial size (no?)?
#

	##### Germination
	# Create Surv object
	germ <- with(subset(data, !is.na(data$MinGermDay)), Surv(MinGermDay, MaxGermDay,	
		type = "interval2"))

	# lognormal produced best fit to data (ID = Population effect)
	fit1 <- survreg(germ ~ Population + frailty(Family, sparse = F), 
		data = subset(data, !is.na(data$MinGermDay)), dist = "lognormal")
	anova(fit1) # significant family and population effects
	
	# Add population germination effects to climate data.frame for subsequent analysis
	focClim$DoG <- c(coef(fit1)[1], coef(fit1)[1] + coef(fit1)[2:16])
	
	# Figure of germination time by latitude
	pdf("ms/Figures/Figure_DoG_Lat.pdf", 5, 5)
	par(mar = c(5, 5, 1, 1))
	plot(focClim$Latitude, c(coef(fit1)[1], coef(fit1)[1] + coef(fit1)[2:16]), 
		xlim = c(30, 45), ylim = c(1.65, 2.65), xlab = "Latitude of Origin", 
		ylab = "Day of Germination [log scale]", 
		cex.lab = 1.5, type = "n", axes = F, frame.plot = T)
	lowCI <- c(confint(fit1)[1, 1], coef(fit1)[1] + confint(fit1)[2:16, 1])
	uppCI <- c(confint(fit1)[1, 2], coef(fit1)[1] + confint(fit1)[2:16, 2])
	arrows(focClim$Latitude, lowCI, focClim$Lat, uppCI, length = 0.05, angle = 90, code = 3)
	points(focClim$Latitude, c(coef(fit1)[1], coef(fit1)[1] + coef(fit1)[2:16]), 
		col = "white", bg = "black", pch = 21, cex = 1.5)
	
	axis(1, at = seq(30, 45, 5), labels = expression(30*degree, 35*degree, 40*degree,
		45*degree), lwd = 0, lwd.ticks = 1)	
	axis(2, lwd = 0, lwd.ticks = 1, las = 1, at = log(seq(6, 14, 2)),
		labels = seq(6, 14, 2))	
	
	fit <- lm(c(coef(fit1)[1], coef(fit1)[1] + coef(fit1)[2:16]) ~ Latitude, data = focClim)
	points(range(focClim$Latitude), predict(fit, new = data.frame(Latitude = 
		range(focClim$Latitude))), type = "l", lwd = 2)
	
	text(grconvertX(0.05, "npc"), grconvertY(0.95, "npc"), adj = c(0, 1),
		labels = bquote(italic(P) == .(round(summary(fit)[[4]][2, 4], 3))),)
	
	dev.off()

	##### Initial size (LLL on May 12)
	data$low <- ifelse(is.na(data$LLL_0512), -Inf, log(data$LLL_0512))
	data$high <- ifelse(is.na(data$LLL_0512), log(0.24), log(data$LLL_0512))
	set.seed(1121014)
	mm1 <- MCMCglmm(cbind(low, high) ~ Population * WaterTrt * TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm2 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt + TempTrt + 
		Population:WaterTrt + Population:TempTrt + WaterTrt:TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm2$DIC - mm1$DIC # drop 3-way interaction
	
	# Test for WaterTrt:TempTrt
	mm3 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt + TempTrt + 
		Population:WaterTrt + Population:TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm3$DIC - mm2$DIC # WaterTrt:TempTrt n.s.

	# Test for Population:TempTrt
	mm4 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt + TempTrt + 
		Population:WaterTrt + WaterTrt:TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm4$DIC - mm2$DIC # Population:TempTrt n.s.

	# Test for Population:WaterTrt
	mm5 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt + TempTrt + 
		Population:TempTrt + WaterTrt:TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm5$DIC - mm2$DIC # Population:WaterTrt n.s.

	##### Drop Population:TempTrt and retest (keep mm4)
	# Test for WaterTrt:TempTrt
	mm6 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt + TempTrt + 
		Population:WaterTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm6$DIC - mm4$DIC # n.s.
	
	# Test for Population:WaterTrt
	mm7 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt + TempTrt + 
		WaterTrt:TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm7$DIC - mm4$DIC # n.s.

	##### Drop Population:WaterTrt and retest (keep mm7)
	# Test for WaterTrt:TempTrt
	mm8 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt + TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm8$DIC - mm7$DIC # n.s.

	##### Drop WaterTrt:TempTrt and retest (keep mm8)
	# Test for TempTrt
	mm9 <- MCMCglmm(cbind(low, high) ~ Population + WaterTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm9$DIC - mm8$DIC # n.s.

	# Test for WaterTrt
	mm10 <- MCMCglmm(cbind(low, high) ~ Population + TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm10$DIC - mm8$DIC # n.s.

	# Test for Population
	mm11 <- MCMCglmm(cbind(low, high) ~ WaterTrt + TempTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm11$DIC - mm8$DIC # n.s.

	##### Drop TempTrt and retest (keep mm9)
	# Test for WaterTrt
	mm12 <- MCMCglmm(cbind(low, high) ~ Population, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm12$DIC - mm9$DIC # n.s.

	# Test for Population
	mm13 <- MCMCglmm(cbind(low, high) ~ WaterTrt, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm13$DIC - mm9$DIC # n.s.

	##### Drop WaterTrt and retest (keep mm12)
	# Test for Population
	mm14 <- MCMCglmm(cbind(low, high) ~ 1, 
		random = ~ Family, data = data, family = "cengaussian", nitt = 1.1e5, thin = 1e2,
		burnin = 1e4)
	mm14$DIC - mm12$DIC # n.s.

	mm15 <- MCMCglmm(cbind(low, high) ~ 1, data = data, family = "cengaussian", 
		nitt = 1.1e5, thin = 1e2, burnin = 1e4)
	mm15$DIC - mm14$DIC # n.s.

	#mm14 is best (almost lowest DIC)
	# InitSizeModels <- list(mm1, mm2, mm3, mm4, mm5, mm6, mm7, mm8, mm9, mm10, mm11,
		# mm12, mm13, mm14, mm15)
	# save(InitSizeModels, file = "Analysis/InitSizeModels.RData")
	load(file = "Analysis/InitSizeModels.RData")
	# Make supp table of results?
	
	rownames <- c("Population + Water + Temperature + Population:Water + Population:Temperature + Water:Temperature + Population:Water:Temperature",
	"Population + Water + Temperature + Population:Water + Population:Temperature + Water:Temperature",
	"Population + Water + Temperature + Population:Water + Population:Temperature",
	"Population + Water + Temperature + Population:Water + Water:Temperature",
	"Population + Water + Temperature + Population:Temperature + Water:Temperature",
	"Population + Water + Temperature + Population:Water",
	"Population + Water + Temperature + Water:Temperature",
	"Population + Water + Temperature", "Population + Water",
	"Population + Temperature", "Water + Temperature", "Population", "Water", "-", "-")

	colnames <- c("Random", "DIC")
	dat <- data.frame(Random = c(rep("Family", 14), "-"), 
		DIC = round(sapply(InitSizeModels, function(X) X$DIC), 1),
		stringsAsFactors = F)

	exportTable(file = "ms/Tables/Table_InitialSize.txt", data = dat, 
		colnames = colnames, rownames = rownames, h1 = "Model")

### Question 1: Variation in 'intrinsic' traits (significant main effect of Population)
### Question 2: Variation in plasticity (significant Treatment x Population interaction)

	# Leaf expansion rate (= LLL Growth rate)
	# Using Absolute growth rates for now as this seems most sensible. Need to test that this correlates with final biomass size.
	
	# Mixed model ANOVA with family as random factor fit using step down procedure
	mm <- lmer(lll_AbsGrowth ~ AvgGermDay + ID * TempTrt * WaterTrt + (1|Family),
		data = data)
	fit1 <- step(mm, reduce.random = F) #includes pop, temp, h2o, population x water
	# fit1 <- step(mm, ddf = "Kenward-Roger", reduce.random = F) # same result using Satterthwaite ddf

	# Dyanamic table of mm output for ms.tex
	tab <- fit1$anova.table[c("AvgGermDay", "Population", "TempTrt", "WaterTrt", 
		"Population:TempTrt", "Population:WaterTrt", "TempTrt:WaterTrt", 
		"Population:TempTrt:WaterTrt"), ]
	
	# Modify P-values for scientific notation
	p <- -expm1(pf(tab$F, tab$NumDF, tab$DenDF, log = T))
	oom <- ceiling(-log10(p)) # order of magnitude
	p <- mapply(function(X1, X2) if(X1 > 0.01){round(X1, 2)} else {
		paste(round(X1 * 10 ^ X2, 2), sprintf("$\\times10^{-%s}$", X2))}, p, oom)
	
	dat <- data.frame(SS = round(tab[, "Sum Sq"], 1), MS = round(tab[, "Mean Sq"], 1),
		df1 = tab[, "NumDF"], df2 = round(tab[, "DenDF"], 1), 
		F = round(tab[, "F.value"], 1), P = p,
		stringsAsFactors = F)
	colnames <- c("SS", "MS", "df1", "df2", "\\em{F}", "$P$")
	rownames <- c("Day of Germination", "Population", "Temperature", "Water", 
		"Population $\\times$ Temperature", "Population $\\times$ Water", 
		"Temperature $\\times$ Water", 
		"Population $\\times$ Temperature $\\times$ Water")
	exportTable(file = "ms/Tables/Table_lllGrowthAnova.txt", data = dat, 
		colnames = colnames, rownames = rownames)
			
	# Plot latitudinal variation in Intercept
	
	# best fit model based on stepdown elimination procedure
	fit1 <- lmer(lll_AbsGrowth ~ AvgGermDay + ID + TempTrt + WaterTrt + ID:WaterTrt + 
		(1|Family), data = data)
		
	# Intercept: least-square coefficients ('betas') and 95% CIs
	betas <- lsmeans(fit1)
	
	# Rows with Population coefficients
	x <- grep("ID  [a-zA-Z]+", rownames(betas$lsmeans.table))

	# Add population leaf elongation rate to climate data.frame for subsequent analysis
	focClim$LER <- betas$lsmeans.table$Estimate[x]

	pdf("ms/Figures/Figure_LLL_Lat.pdf", 5, 5)
	par(mar = c(5, 5, 1, 1))
	plot(focClim$Latitude, betas$lsmeans.table$Estimate[x], xlim = c(30, 45), 
		ylim = c(1.5, 3.1), xlab = "Latitude of Origin", 
		ylab = expression(paste("Leaf elongation rate (", mm~day^{-1}, ")")), 
		las = 1, cex.lab = 1.5, type = "n", axes = F, frame.plot = T)
	arrows(focClim$Lat, betas$lsmeans.table[x, "Lower CI"], focClim$Lat, 
		betas$lsmeans.table[x, "Upper CI"], angle = 90, length = 0.05, code = 3)
	points(focClim$Lat, betas$lsmeans.table$Estimate[x], col = "white", bg = "black",
		pch = 21, cex = 1.5)
	
	axis(1, at = seq(30, 45, 5), labels = expression(30*degree, 35*degree, 40*degree,
		45*degree), lwd = 0, lwd.ticks = 1)	
	axis(2, lwd = 0, lwd.ticks = 1, las = 1)	
	
	fit <- lm(betas$lsmeans.table$Estimate[x] ~ Latitude, data = focClim)

	points(range(focClim$Latitude), predict(fit, new = data.frame(Latitude = 
		range(focClim$Latitude))), type = "l", lwd = 2)

	text(grconvertX(0.95, "npc"), grconvertY(0.95, "npc"), adj = c(1, 1),
		labels = bquote(italic(P) == .(round(summary(fit)[[4]][2, 4], 3))))
	
	dev.off()
			
	# Stem elongation rate (= Height Growth rate)

	# Mixed model ANOVA with family as random factor fit using step down procedure
	mm <- lmer(height_AbsGrowth ~ AvgGermDay + Population * TempTrt * WaterTrt + 
		(1|Family), data = data)
	fitSER <- step(mm, reduce.random = F) #includes pop, temp, h2o, temp x water
	# fit2 <- step(mm, ddf = "Kenward-Roger", reduce.random = F) # same result using Satterthwaite ddf

	# Dyanamic table of mm output for ms.tex
	tab <- fit1$anova.table[c("AvgGermDay", "Population", "TempTrt", "WaterTrt",
		"Population:TempTrt", "Population:WaterTrt", "TempTrt:WaterTrt", 
		"Population:TempTrt:WaterTrt"), ]
	
	# Modify P-values for scientific notation
	p <- -expm1(pf(tab$F, tab$NumDF, tab$DenDF, log = T))
	oom <- ceiling(-log10(p)) # order of magnitude
	p <- mapply(function(X1, X2) if(X1 > 0.01){round(X1, 2)} else {
		paste(round(X1 * 10 ^ X2, 2), sprintf("$\\times10^{-%s}$", X2))}, p, oom)
	
	dat <- data.frame(SS = round(tab[, "Sum Sq"], 1), MS = round(tab[, "Mean Sq"], 1),
		df1 = tab[, "NumDF"], df2 = round(tab[, "DenDF"], 1), 
		F = round(tab[, "F.value"], 1), P = p,
		stringsAsFactors = F)
	colnames <- c("SS", "MS", "df1", "df2", "\\em{F}", "$P$")
	rownames <- c("Day of Germination", "Population", "Temperature", "Water", 
		"Population $\\times$ Temperature", "Population $\\times$ Water", 
		"Temperature $\\times$ Water", 
		"Population $\\times$ Temperature $\\times$ Water")
	exportTable(file = "ms/Tables/Table_heightGrowthAnova.txt", data = dat, 
		colnames = colnames, rownames = rownames)

	# Plot latitudinal variation in Intercept and Population:WaterTrt
	
	# best fit model based on stepdown elimination procedure
	fit1 <- lmer(height_AbsGrowth ~ AvgGermDay + ID + TempTrt + WaterTrt + 
		TempTrt:WaterTrt + (1|Family), data = data)

	# Intercept: least-square coefficients ('betas') and 95% CIs
	betas <- lsmeans(fit1)
	
	# Rows with Population coefficients
	x <- grep("ID  [a-zA-Z]+", rownames(betas$lsmeans.table))

	# Add population stem elongation rate to climate data.frame for subsequent analysis
	focClim$SER <- betas$lsmeans.table$Estimate[x]

	pdf("ms/Figures/Figure_Height_Lat.pdf", 5, 5)
	par(mar = c(5, 5, 1, 1))
	plot(focClim$Lat, betas$lsmeans.table$Estimate[x], xlim = c(30, 45), 
		ylim = c(0.75, 1.65), xlab = "Latitude of Origin", 
		ylab = expression(paste("Stem elongation rate (", mm~day^{-1}, ")")), 
		cex.lab = 1.5, type = "n", axes = F, frame.plot = T)
	arrows(focClim$Lat, betas$lsmeans.table[x, "Lower CI"], focClim$Lat, 
		betas$lsmeans.table[x, "Upper CI"], angle = 90, length = 0.05, code = 3)
	points(focClim$Lat, betas$lsmeans.table$Estimate[x], col = "white", bg = "black",
		pch = 21, cex = 1.5)
	
	axis(1, at = seq(30, 45, 5), labels = expression(30*degree, 35*degree, 40*degree,
		45*degree), lwd = 0, lwd.ticks = 1)	
	axis(2, lwd = 0, lwd.ticks = 1, las = 1)	
	
	fit <- lm(betas$lsmeans.table$Estimate[x] ~ Lat, data = focClim)
	points(range(focClim$Lat), predict(fit, new = data.frame(Lat = range(focClim$Lat))),
		type = "l", lwd = 2)
	
	text(grconvertX(0.95, "npc"), grconvertY(0.95, "npc"), adj = c(1, 1),
		labels = bquote(italic(P) == .(round(summary(fit)[[4]][2, 4], 3))),)
	
	dev.off()
			
  # Photosynthetic rate
  mm <- lmer(Photo ~ Population * TempTrt * WaterTrt + (1|Family),
    data = subset(data, !is.na(data$Photo)))
	fitPhoto <- step(mm, reduce.random = F)
  
	mm <- lmer(resPhoto ~ Population * TempTrt * WaterTrt + (1|Family),
	           data = subset(data, !is.na(data$Photo))[-52, ])
	fitIntrnPhoto <- step(mm, reduce.random = F)
  
	# Morality
	fitMort <- glmer(DiedFromStress ~ Population * TempTrt * WaterTrt + (1|Family) + (1|indiv),
    data = subset(data, !is.na(data$DiedFromStress)), family = "binomial")

  fitMort <- glm(DiedFromStress ~ Population * TempTrt * WaterTrt,
	                 data = subset(data, !is.na(data$DiedFromStress)), family = "quasibinomial")

  Anova(fitMort, type = 2)
  
  
  #### TO DO: explain variation in plasticity...
# PART II: SELECTIVE AGENT

	# Load vector of exploratory climate variables used in climate to latitude model search
	load("~/Google Drive/CardLocalAdaptation/Data/Climate/Robjects/Exploratory_climate_variables.RData") #climVars
	
	# Run VSURF
	set.seed(022918)
	linmod.DoG <- as.formula(paste("DoG ~ ", paste(climVars, collapse = " + "), sep = ""))
	linmod.LER <- as.formula(paste("LER ~ ", paste(climVars, collapse = " + "), sep = ""))
	linmod.SER <- as.formula(paste("SER ~ ", paste(climVars, collapse = " + "), sep = ""))
	fitDoG <- VSURF.parallel(linmod.DoG, data = focClim[, c("DoG", climVars)], ncores = 3)
	fitLER <- VSURF.parallel(linmod.LER, data = focClim[, c("LER", climVars)], ncores = 3)
	fitSER <- VSURF.parallel(linmod.SER, data = focClim[, c("SER", climVars)], ncores = 3)
	climVars[fitDoG$varselect.pred]
	climVars[fitLER$varselect.pred]
	climVars[fitSER$varselect.pred]
	
	# save(fitDoG, file = "Analysis/clim2DoG_fit.RData")
	# save(fitLER, file = "Analysis/clim2LER_fit.RData")
	# save(fitSER, file = "Analysis/clim2SER_fit.RData")
	load(file = "Analysis/clim2DoG_fit.RData")
	load(file = "Analysis/clim2LER_fit.RData")
	load(file = "Analysis/clim2SER_fit.RData")
	
	# Figure for comparison climate 2 latitude analysis
	
	load(file = paste(path.obj, "/clim2lat_fit.RData", sep = "")) # fit
	load(file = paste(path.obj, "/clim2lat_fit_t1.RData", sep = "")) # fit_t1
	
	# Important climatic variables for predicting latitude
	impVars <- climVars[fit_t1$varselect.pred]
	nImpLat <- length(fit_t1$varselect.pred)
	nImpDoG <- length(fitDoG$varselect.pred)
	nImpLER <- length(fitLER$varselect.pred)
	nImpSER <- length(fitSER$varselect.pred)
	
	# Figure of variable importance
	pdf("~/Google Drive/CardLocalAdaptation/ms/Figures/Figure_clim2traits_imp.pdf", 4, 5)
	par(mar = c(5, 1, 1, 6))
	
	# Day of Germination
	plot(0, 0, type = "n", xlim = c(0, fitDoG$ord.imp$x[1] * 1000), ylim = c(1, nImpLat), 
		axes = F, frame.plot = T, xlab = "Importance", ylab = "", cex.lab = 1.5)
	for (i in (nImpLat:(nImpLat - nImpDoG + 1))) lines(x = c(-100, 100), y = c(i, i), 
		col = "grey")
	points(1000 * fitDoG$ord.imp$x[1:nImpDoG], nImpLat:(nImpLat - nImpDoG + 1), pch = 21, 
		col = "black", bg = "white", cex = 1.5)
	mtext(climVars[fitDoG$varselect.pred], side = 4, at = nImpLat:(nImpLat - nImpDoG + 1), 
		las = 1, line = 1, font = ifelse(fitDoG$varselect.pred %in% fit_t1$varselect.pred,
		2, 1))
	axis(1, lwd = 0, lwd.ticks = 1, labels = seq(0, 2, 0.5) / 1e3, at = seq(0, 2, 0.5))

	# Leaf expansion rate
	plot(0, 0, type = "n", xlim = c(0, fitLER$ord.imp$x[1] * 1000), ylim = c(1, nImpLat), 
		axes = F, frame.plot = T, xlab = "Importance", ylab = "", cex.lab = 1.5)
	for (i in (nImpLat:(nImpLat - nImpLER + 1))) lines(x = c(-100, 100), y = c(i, i), 
		col = "grey")
	points(1000 * fitLER$ord.imp$x[1:nImpLER], nImpLat:(nImpLat - nImpLER + 1), pch = 21, 
		col = "black", bg = "white", cex = 1.5)
	mtext(climVars[fitLER$varselect.pred], side = 4, at = nImpLat:(nImpLat - nImpLER + 1), 
		las = 1, line = 1, font = ifelse(fitLER$varselect.pred %in% fit_t1$varselect.pred,
		2, 1))
	axis(1, lwd = 0, lwd.ticks = 1, labels = seq(0, 1.5, 0.5) / 1e3, at = seq(0, 1.5, 0.5))

	# Stem elongation rate
	plot(0, 0, type = "n", xlim = c(0, fitSER$ord.imp$x[1] * 1000), ylim = c(1, nImpLat), 
		axes = F, frame.plot = T, xlab = "Importance", ylab = "", cex.lab = 1.5)
	for (i in (nImpLat:(nImpLat - nImpSER + 1))) lines(x = c(-100, 100), y = c(i, i), 
		col = "grey")
	points(1000 * fitSER$ord.imp$x[1:nImpSER], nImpLat:(nImpLat - nImpSER + 1), pch = 21, 
		col = "black", bg = "white", cex = 1.5)
	mtext(climVars[fitSER$varselect.pred], side = 4, at = nImpLat:(nImpLat - nImpSER + 1), 
		las = 1, line = 1, font = ifelse(fitSER$varselect.pred %in% fit_t1$varselect.pred,
		2, 1))
	axis(1, lwd = 0, lwd.ticks = 1, labels = seq(0, 3, 0.5) / 1e3, at = seq(0, 3, 0.5))

	dev.off()
	
	
#### code for multimodal avg if I decide to use this:

library("MuMIn")

# test all two parameter combinations for germination (DoG), leaf elongation rate (LER), and stem elongation rate (SER), and do multimodel average
focClim$DoG <- c(coef(fit1)[1], coef(fit1)[1] + coef(fit1)[2:16])
focClim$LER <- betas$lsmeans.table$Estimate[x]
focClim$SER <- betas$lsmeans.table$Estimate[x]

# should latitude be included in predictors???
predCombos <- combn(c(2, 6:57), 2)
par2Formulae <- apply(predCombos, 2, function(X)
{
	as.formula(paste("DoG ~ ", paste(colnames(focClim)[X], collapse = " + "), sep = ""))
})

fitDoG <- lapply(par2Formulae, lm, data = focClim)

bestModel <- which.min(sapply(fitDoG, AIC))
par2Formulae[[bestModel]]
summary(fitDoG[[bestModel]])

# models with delta(AIC) < 2
dAIC <- sapply(fitDoG, AIC) - min(sapply(fitDoG, AIC))
par2Formulae[[which(dAIC <= 2)[1]]]


#####

par2Formulae <- apply(predCombos, 2, function(X)
{
	as.formula(paste("LER ~ ", paste(colnames(focClim)[X], collapse = " + "), sep = ""))
})

fitLER <- lapply(par2Formulae, lm, data = focClim)

bestModel <- which.min(sapply(fitLER, AIC))
par2Formulae[[bestModel]]
summary(fitLER[[bestModel]])

# models with delta(AIC) < 2
dAIC <- sapply(fitLER, AIC) - min(sapply(fitLER, AIC))
par2Formulae[[which(dAIC <= 2)[3]]]

# multimodel inference
mmAvg <- model.sel(fitLER)
importance(mmAvg)


#####

par2Formulae <- apply(predCombos, 2, function(X)
{
	as.formula(paste("SER ~ ", paste(colnames(focClim)[X], collapse = " + "), sep = ""))
})

fitSER <- lapply(par2Formulae, lm, data = focClim)

bestModel <- which.min(sapply(fitSER, AIC))
par2Formulae[[bestModel]]
summary(fitSER[[bestModel]])

# models with delta(AIC) < 2
dAIC <- sapply(fitSER, AIC) - min(sapply(fitSER, AIC))
par2Formulae[[which(dAIC <= 2)[1]]]

# multimodel inference
mmAvg <- model.avg(fitSER)
importance(mmAvg)
