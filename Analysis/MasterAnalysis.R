#
#	Master Analysis of M. cardinalis Growth Chamber Experiment (Spring 2014)
#
#	Title:	 Local Physiological adaptaiton in Mimulus cardinalist
#	Authors: Christopher D. Muir and Amy L. Angert
#	Journal: Evolution?
#

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
	
#
##	Set working directory
#

	setwd("~/Google Drive/CardLocalAdaptation")

#
##	Source custom functions
#

	source("Analysis/functions.R")
	
#
##	Population Environmental data
#
	# This needs to be updated with PRISM, other data, etc. Lat, lon should be correct though
	popenv <- read.csv("Data/PopulationEnvData.csv")
	
#
##	Read in processed datasets
#

	data <- read.csv("Data/MasterDatasheet_out.csv")
	
### JUST USING GROWTH NOW UNTIL OTHER DATA ARE READY
### Analysis 1: Variation in 'intrinsic' traits (significant main effect of Population)
### Analysis 2: Variation in plasticity (significant Treatment x Population interaction)

	# LLL Growth rate
	# Using Absolute growth rates for now as this seems most sensible. Need to test that this correlates with final biomass size.
	
	# approach 1: stepAIC, no Family effect
	fit1 <- stepAIC(lm(lll_AbsGrowth ~ Population * TempTrt * WaterTrt, data = data))
	Anova(fit1, type = 2)

	# approach 2: mixed model with family	
	mm <- lmer(lll_AbsGrowth ~ Population * TempTrt * WaterTrt + (1|Family),
		data = data)
	fit1 <- step(mm) #includes pop, temp, h2o, population x water
	# fit1 <- step(mm, ddf = "Kenward-Roger") # same result using Satterthwaite ddf

	# Dyanamic table of mm output for ms.tex
	tab <- fit1$anova.table[c("Population", "TempTrt", "WaterTrt", "Population:TempTrt",
		"Population:WaterTrt", "TempTrt:WaterTrt", "Population:TempTrt:WaterTrt"), ]
	
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
	rownames <- c("Population", "Temperature", "Water", 
		"Population $\\times$ Temperature", "Population $\\times$ Water", 
		"Temperature $\\times$ Water", 
		"Population $\\times$ Temperature $\\times$ Water")
	exportTable(file = "ms/Tables/Table_lllGrowthAnova.txt", data = dat, 
		colnames = colnames, rownames = rownames)
			
	# Plot latitudinal variation
	
	#### NEED TO DO ####

	# Height Growth rate
	# Using Absolute growth rates for now as this seems most sensible. Need to test that this correlates with final biomass size.

	# approach 1: stepAIC, no Family effect
	fit1 <- stepAIC(lm(height_AbsGrowth ~ Population * TempTrt * WaterTrt, data = data))
	Anova(fit1, type = 2)

	# approach 2: mixed model with family	
	mm <- lmer(height_AbsGrowth ~ Population * TempTrt * WaterTrt + (1|Family),
		data = data)
	fit1 <- step(mm, reduce.random = F) #includes pop, temp, h2o, temp x water
	# fit1 <- step(mm, ddf = "Kenward-Roger", reduce.random = F) # same result using Satterthwaite ddf

	# Dyanamic table of mm output for ms.tex
	tab <- fit1$anova.table[c("Population", "TempTrt", "WaterTrt", "Population:TempTrt",
		"Population:WaterTrt", "TempTrt:WaterTrt", "Population:TempTrt:WaterTrt"), ]
	
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
	rownames <- c("Population", "Temperature", "Water", 
		"Population $\\times$ Temperature", "Population $\\times$ Water", 
		"Temperature $\\times$ Water", 
		"Population $\\times$ Temperature $\\times$ Water")
	exportTable(file = "ms/Tables/Table_heightGrowthAnova.txt", data = dat, 
		colnames = colnames, rownames = rownames)






#
# DFA/PCA
#

	fit <- lda(Population ~ b0_lll + b1_height,
		data = subset(data, data$TempTrt == "Hot"), CV = T)
	ct <- table(data$Population[data$TempTrt == "Hot"], fit$class)
	diag(prop.table(ct, 1))
	# total percent correct
	sum(diag(prop.table(ct)))

# scratch

tapply(lll.data$LLL, lll.data$indiv, max)
plot(ranef(mm)$indiv[, 1], tapply(lll.data$LLL, lll.data$indiv, max))
plot(ranef(mm)$indiv[, 3], tapply(height.data$height, height.data $indiv, max))