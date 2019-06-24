
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
# load the dataset
reaction <- read.table("ReactionNoStapel.csv", header = TRUE, sep = ",")
covariance <- read.table("CovarianceReaction.csv", header = FALSE, sep = ",")
#for some reason covariance is read with an extra line at the end, remove it
V <- data.frame(covariance)
V <- V[-c(186),]
V<-data.matrix(V)
V <- forceSymmetric(V)
reaction$WT <- 1/reaction$InvWt_UpDn

# Set the moderator name
reaction$moderator <- reaction$Task
# Multivariate approach
reaction$Task <- relevel(factor(reaction$Task), ref="1")
reaction$DV_Code <- relevel(factor(reaction$DV_Code),ref="3")
moderator <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code*moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
addmargins(table(reaction$DV_Code, reaction$moderator))
moderator
# may need some changes here
anova(moderator, btt=2:5)
anova(moderator, btt=7:10)
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(moderator, sigma2=1)
profile(moderator, sigma2=2)

# Then to get the predicted levels for each part, do this...
noveltypred <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code:Known.dimension- 1, random=~1 | Study_Code, data=reaction, subset = (Task==1))
noveltypred

# Funnel plot for subgroups
noveltyKnown <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code, random = ~ 1 | Study_Code, data=reaction, subset = (Known.dimension==1 & DV_Code ==3))

noveltyKnown <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code, random = ~ 1 | Study_Code, data=reaction, subset = (Known.dimension==1))
novelty <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code*Known.dimension, random = ~ 1 | Study_Code, data=reaction, subset = (Task==1))
funnel(overalluni)
overallunireg <- regtest.rma(overalluni)
print.regtest.rma(overallunireg)

# trim-and-fill (will have one for each level)
known_ability <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==3) & (Known.dimension ==1))
trimfill(known_ability, "right")
funnel(known_ability)
knownabilityregtest <- regtest.rma(known_ability)
print.regtest.rma(knownabilityregtest)

novel_ability <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==3) & (Known.dimension ==2))
trimfill(novel_ability, "right")
funnel(novel_ability)
novelabilityregtest <- regtest.rma(novel_ability)
print.regtest.rma(novelabilityregtest)

known_affect <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==1) & (Known.dimension ==1))
trimfill(known_affect, "right")
funnel(known_affect)
knownaffectregtest <- regtest.rma(known_affect)
print.regtest.rma(knownaffectregtest)

novel_affect <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==1) & (Known.dimension ==2))
trimfill(novel_affect, "right")
funnel(novel_affect)
novelaffectregtest <- regtest.rma(novel_affect)
print.regtest.rma(novelaffectregtest)

known_se <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==2) & (Known.dimension ==1))
trimfill(known_se, "right")
funnel(known_se)
knownseregtest <- regtest.rma(known_se)
print.regtest.rma(knownseregtest)

known_beh <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==4) & (Known.dimension ==1))
trimfill(known_beh, "right")
funnel(known_beh)
knownbehregtest <- regtest.rma(known_beh)
print.regtest.rma(knownbehregtest)

novel_beh <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==4) & (Known.dimension ==2))
trimfill(novel_beh, "right")
funnel(novel_beh)
novelbehregtest <- regtest.rma(novel_beh)
print.regtest.rma(novelbehregtest)

known_psat <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==5) & (Known.dimension ==1))
trimfill(known_psat, "right")
funnel(known_psat)
knownpsatregtest <- regtest.rma(known_psat)
print.regtest.rma(knownpsatregtest)

novel_psat <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==5) & (Known.dimension ==2))
trimfill(novel_psat, "right")
funnel(novel_psat)
novelpsatregtest <- regtest.rma(novel_psat)
print.regtest.rma(novelpsatregtest)


# Robust version
library(clubSandwich)
novelty_robust <- coef_test(novelty, vcov = "CR2")
novelty_robust

# Vevea & Woods
library(weightr)
abilityWoods <- reaction[reaction$DV_Code == 3,]
abknownWoods <- abilityWoods[abilityWoods$Known.dimension==1,]
abnovelWoods <- abilityWoods[abilityWoods$Known.dimension==2,]
abknVW2tailmod <-weightfunct(abknownWoods$ES_Up_Down, abknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abknVW2tailmod[[2]]$par
abknVW2tailsevere <-weightfunct(abknownWoods$ES_Up_Down, abknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abknVW2tailsevere[[2]]$par
abnoVW2tailmod <-weightfunct(abnovelWoods$ES_Up_Down, abnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abnoVW2tailmod[[2]]$par
abnoVW2tailsevere <-weightfunct(abnovelWoods$ES_Up_Down, abnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abnoVW2tailsevere[[2]]$par

afWoods <- reaction[reaction$DV_Code == 1,]
afknownWoods <- afWoods[afWoods$Known.dimension==1,]
afnovelWoods <- afWoods[afWoods$Known.dimension==2,]
afknVW2tailmod <-weightfunct(afknownWoods$ES_Up_Down, afknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
afknVW2tailmod[[2]]$par
afknVW2tailsevere <-weightfunct(afknownWoods$ES_Up_Down, afknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
afknVW2tailsevere[[2]]$par
afnoVW2tailmod <-weightfunct(afnovelWoods$ES_Up_Down, afnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
afnoVW2tailmod[[2]]$par
afnoVW2tailsevere <-weightfunct(afnovelWoods$ES_Up_Down, afnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
afnoVW2tailsevere[[2]]$par

seWoods <- reaction[reaction$DV_Code == 2,]
seknownWoods <- seWoods[seWoods$Known.dimension==1,]
senovelWoods <- seWoods[seWoods$Known.dimension==2,]
seknVW2tailmod <-weightfunct(seknownWoods$ES_Up_Down, seknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
seknVW2tailmod[[2]]$par
seknVW2tailsevere <-weightfunct(seknownWoods$ES_Up_Down, seknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
seknVW2tailsevere[[2]]$par
senoVW2tailmod <-weightfunct(senovelWoods$ES_Up_Down, senovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
senoVW2tailmod[[2]]$par
senoVW2tailsevere <-weightfunct(senovelWoods$ES_Up_Down, senovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
senoVW2tailsevere[[2]]$par

beWoods <- reaction[reaction$DV_Code == 4,]
beknownWoods <- beWoods[beWoods$Known.dimension==1,]
benovelWoods <- beWoods[beWoods$Known.dimension==2,]
beknVW2tailmod <-weightfunct(beknownWoods$ES_Up_Down, beknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
beknVW2tailmod[[2]]$par
beknVW2tailsevere <-weightfunct(beknownWoods$ES_Up_Down, beknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
beknVW2tailsevere[[2]]$par
benoVW2tailmod <-weightfunct(benovelWoods$ES_Up_Down, benovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
benoVW2tailmod[[2]]$par
benoVW2tailsevere <-weightfunct(benovelWoods$ES_Up_Down, benovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
benoVW2tailsevere[[2]]$par

psWoods <- reaction[reaction$DV_Code == 5,]
psknownWoods <- psWoods[psWoods$Known.dimension==1,]
psnovelWoods <- psWoods[psWoods$Known.dimension==2,]
psknVW2tailmod <-weightfunct(psknownWoods$ES_Up_Down, psknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
psknVW2tailmod[[2]]$par
psknVW2tailsevere <-weightfunct(psknownWoods$ES_Up_Down, psknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
psknVW2tailsevere[[2]]$par
psnoVW2tailmod <-weightfunct(psnovelWoods$ES_Up_Down, psnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
psnoVW2tailmod[[2]]$par
psnoVW2tailsevere <-weightfunct(psnovelWoods$ES_Up_Down, psnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
psnoVW2tailsevere[[2]]$par

#PEESE
knabilitypeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==3) & (Known.dimension==1))
knabilitypeese
noabilitypeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==3) & (Known.dimension==2))
noabilitypeese
knaffectpeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==1) & (Known.dimension==1))
knaffectpeese
noaffectpeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==1) & (Known.dimension==2))
noaffectpeese
knsepeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==2) & (Known.dimension==1))
knsepeese
knbehpeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==4) & (Known.dimension==1))
knbehpeese
nobehpeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==4) & (Known.dimension==2))
nobehpeese
knperfspeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==5) & (Known.dimension==1))
knperfspeese
noperfspeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE", subset = (DV_Code==5) & (Known.dimension==2))
noperfspeese
