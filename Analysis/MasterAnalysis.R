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
	
#
##	Set working directory
#

	setwd("~/Google Drive/CardLocalAdaptation")

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
	fit1 <- lmer(lll_AbsGrowth ~ Population * TempTrt * WaterTrt + (1|Family), data = data)
	
	# 3-way interaction
	fit2 <- update(fit1, . ~ . - Population:TempTrt:WaterTrt)
	KRmodcomp(fit1, fit2) # n.s.

	# Population x Treatment interactions
	fit3 <- update(fit2, . ~ . - Population:TempTrt)
	fit4 <- update(fit2, . ~ . - Population:WaterTrt)
	KRmodcomp(fit2, fit3) # n.s.
	KRmodcomp(fit2, fit4) # marginally significant. keep.
	
	# Main effects
	fit5 <- update(fit4, . ~ . - Population)
	fit6 <- update(fit4, . ~ . - TempTrt)
	fit7 <- update(fit4, . ~ . - WaterTrt)
	KRmodcomp(fit4, fit5) # n.s.
	KRmodcomp(fit4, fit6) # marginally significant. keep.
	KRmodcomp(fit4, fit7) # marginally significant. keep.

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