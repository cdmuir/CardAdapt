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
	# fit1 <- step(fit1, ddf = "Kenward-Roger") # same result using Satterthwaite ddf

	# Dyanamic table of mm output for ms.tex
	fit1$anova.table
	exportTable(file = "ms/Tables/Table_lllGrowthAnova.txt", data = data.frame(rnorm(5), rnorm(5)), colnames = c("Test1", "Test2"), rownames = LETTERS[1:5])


	# modify code to print anova tables in latex
	mapply(function(X, Y)
	{
		tmp <- anova(X, fit1a)
		p <- tmp[, "Pr(>F)"][2]
		oom <- ceiling(-log10(p)) # order of magnitude
		p <- if(p > 0.001){round(p, 2)}else{paste(round(p * 10 ^ oom, 2), 
			sprintf("$\\times10^{-%s}$", oom))}
		cat(sprintf('%s & %s & %s & %s & %s & %s \\\\', Y, tmp$Df[2], 
			round(tmp[, "Sum of Sq"][2], 3), round(tmp[, "Sum of Sq"][2] / tmp$Df[2], 3), 
			round(tmp$F[2], 2), p), file = table2)
	}, listFita, predVar)
		
	
	# Plot latitudinal variation
	

	# Height Growth rate
	# Using Absolute growth rates for now as this seems most sensible. Need to test that this correlates with final biomass size.

		
	fit1 <- stepAIC(lm(height_AbsGrowth ~ Population * TempTrt * WaterTrt, data = data))
	Anova(fit1, type = 3)
	### Deciding between MANOVA or Multi-response MCMCglmm
	
	# MANOVA [how do I account for family?]
	fit <- manova(cbind(b0_lll, b1_lll, b2_lll) ~ TempTrt * WaterTrt * Population + 
		Error(Family), data[!is.na(data$b0_lll), ])
	summary(fit)
	
	fit <- manova(cbind(b0_lll, b1_lll, b2_lll) ~ Population * TempTrt * WaterTrt,
		data = data[!is.na(data$b0_lll), ])
	Anova(fit, type = 2)

	# Multiresponse model using MCMCglmm [uncertain about prior]
	prior1 <- list(R = list(V = diag(3), nu = 1), G = list(G1 = list(V = diag(2), nu = 2))) 

	mm1 <- MCMCglmm(cbind(b0_lll, b1_lll, b2_lll) ~ TempTrt * WaterTrt, 
		random = ~ us(1 + TempTrt):Family, 
		rcov = ~ cor(trait):units, data = data, family = rep("gaussian", 3), prior = prior1)

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