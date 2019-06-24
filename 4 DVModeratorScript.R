
# Reaction script
##### Preliminary #####################
# clear workspace
# If you need to install the packages install.packages(c("metafor", "clubSandwich", "robumeta", "weightr"))
rm(list=ls())
# store the current directory
initial.dir<-getwd()
# change to the new directory
setwd("/Users/jonathan.gerber/Dropbox/SCMeta_Articles/000 SCMeta Reports/7 Second resubmit/Analyses")
# load the necessary libraries
library(metafor)
# set the output file
#sink(file = "reaction.out", append = TRUE)
#sink()
# load the dataset
reaction <- read.table("ReactionNoStapel.csv", header = TRUE, sep = ",")
covariance <- read.table("CovarianceReaction.csv", header = FALSE, sep = ",")
#for some reason covariance is read with an extra line at the end, remove it
V <- data.frame(covariance)
V <- V[-c(186),]
V<-data.matrix(V)
V <- forceSymmetric(V)
reaction$WT <- 1/reaction$InvWt_UpDn
# Basic univariate type
reaction$DV_Code <- relevel(factor(reaction$DV_Code), ref="3")
uniDV <- rma(ES_Up_Down, WT, mods = ~DV_Code, random = ~1 | Study_Code/Order.for.CMA, data=reaction)
uniDV
predict(uniDV, newmods = rbind(c(0,0,0,0),c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)))

# Multivariate by each DV (multivariate)
reaction$DV_Code <- relevel(factor(reaction$DV_Code), ref="3")
by_DV <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
addmargins(table(reaction$DV_Code))
by_DV
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(by_DV, sigma2=1)
profile(by_DV, sigma2=2)
predict(by_DV, newmods = rbind(c(0,0,0,0),c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)))

# Egger's and trim-and-fill (will have one for each level)
DV_ability <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==3))
trimfill(DV_ability, "right")
funnel(DV_ability)
DVabilityregtest <- regtest.rma(DV_ability)
print.regtest.rma(DVabilityregtest)
ranktest(DV_ability)
hc(DV_ability)

DV_affect <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==1))
ranktest(DV_affect)
funnel(DV_affect)
DVaffectregtest <- regtest.rma(DV_affect)
print.regtest.rma(DVaffectregtest)
trimfill(DV_affect, "right")
hc(DV_affect)

DV_selfesteem <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==2))
ranktest(DV_selfesteem)
funnel(DV_selfesteem)
DVseregtest <- regtest.rma(DV_selfesteem)
print.regtest.rma(DVseregtest)
trimfill(DV_selfesteem, "right")
hc(DV_selfesteem)

DV_behavior <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==4))
ranktest(DV_behavior)
DVbehaviorregtest <- regtest.rma(DV_behavior)
print.regtest.rma(DVbehaviorregtest)
trimfill(DV_behavior, "right")
funnel(DV_behavior)
hc(DV_behavior)

DV_perfsat <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==5))
funnel(DV_perfsat)
ranktest(DV_perfsat)
DVperfsatregtest <- regtest.rma(DV_perfsat)
print.regtest.rma(DVperfsatregtest)
trimfill(DV_perfsat, "left")
hc(DV_perfsat)

# Robust version
library(clubSandwich)
DV_robust <- coef_test(by_DV, vcov = "CR2")
DV_robust

# Vevea & Woods
library(weightr)
abilityWoods <- reaction[reaction$DV_Code == 3,]
overallVW2tailmod <-weightfunct(abilityWoods$ES_Up_Down, abilityWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
overallVW2tailmod[[2]]$par
overallVW2tailsevere <-weightfunct(abilityWoods$ES_Up_Down, abilityWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
overallVW2tailsevere[[2]]$par

DVaffWoods <- reaction[reaction$DV_Code == 1,]
DVaffVW2tailmod <-weightfunct(DVaffWoods$ES_Up_Down, DVaffWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
DVaffVW2tailmod[[2]]$par
DVaffVW2tailsevere <-weightfunct(DVaffWoods$ES_Up_Down, DVaffWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
DVaffVW2tailsevere[[2]]$par

DVseWoods <- reaction[reaction$DV_Code == 2,]
DVseVW2tailmod <-weightfunct(DVseWoods$ES_Up_Down, DVseWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
DVseVW2tailmod[[2]]$par
DVseVW2tailsevere <-weightfunct(DVseWoods$ES_Up_Down, DVseWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
DVseVW2tailsevere[[2]]$par

DVbehWoods <- reaction[reaction$DV_Code == 4,]
DVbehW2tailmod <-weightfunct(DVbehWoods$ES_Up_Down, DVbehWoods$WT, steps=c(0.1, 0.50,.75,1), weights=c(.95,.6,.6,1))
DVbehVW2tailmod[[2]]$par
DVbehVW2tailsevere <-weightfunct(DVbehWoods$ES_Up_Down, DVbehWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
DVbehVW2tailsevere[[2]]$par

perfsatWoods <- reaction[reaction$DV_Code == 5,]
DVperfsatVW2tailmod <-weightfunct(perfsatWoods$ES_Up_Down, perfsatWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
DVperfsatVW2tailmod[[2]]$par
DVperfsatVW2tailsevere <-weightfunct(perfsatWoods$ES_Up_Down, perfsatWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
DVperfsatVW2tailsevere[[2]]$par

DVallVW2tailmod <-weightfunct(abilityWoods$ES_Up_Down, abilityWoods$WT, mods = ~ abilityWoods$DV_Code, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
overallVW2tailmod[[2]]$par

#PEESE
abilitypeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==3))
abilitypeese
confint(abilitypeese)
peeselm <- lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==3))
summary(peeselm)
confint(peeselm)

affectpeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==1))
affectpeese
sepeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==2))
sepeese
behpeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==4))
behpeese
perfspeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==5))
perfspeese

affpeeselm <- lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==1))
summary(affpeeselm)
confint(affpeeselm)
sepeeselm <- lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==2))
summary(sepeeselm)
confint(sepeeselm)
peeselm <- lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==4))
summary(peeselm)
confint(peeselm)
peeselm <- lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==5))
summary(peeselm)
confint(peeselm)
