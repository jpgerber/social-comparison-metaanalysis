
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
reaction$moderator <- reaction$Known.dimension
# Univariate approach
reaction$moderator <- relevel(factor(reaction$moderator), ref="1")
reaction$DV_Code <- relevel(factor(reaction$DV_Code),ref="3")
uni <- rma(ES_Up_Down, WT, mods = ~ DV_Code*moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
uni
predict(uni, newmods = rbind(c(0,0,0,0,0,0,0,0),c(0,0,0,0,1,0,0,0),c(1,0,0,0,0,0,0,0),c(1,0,0,0,1,1,0,0),
                             c(0,1,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0),c(0,0,1,0,1,0,1,0),
                             c(0,0,0,1,0,0,0,0),c(0,0,0,1,1,0,0,1)))

anova(uni, btt=2:5)
anova(uni, btt=7:10)
par(mar=c(1,1,1,1))
par(mfrow=c(1,1))
profile(uni, sigma2=1)

# Multivariate approach
reaction$moderator <- relevel(factor(reaction$moderator), ref="1")
reaction$DV_Code <- relevel(factor(reaction$DV_Code),ref="3")
moderator <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code*moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
addmargins(table(reaction$DV_Code, reaction$moderator))
moderator
predict(moderator, newmods = rbind(c(0,0,0,0,0,0,0,0),c(0,0,0,0,1,0,0,0),c(1,0,0,0,0,0,0,0),c(1,0,0,0,1,1,0,0),
                             c(0,1,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0),c(0,0,1,0,1,0,1,0),
                             c(0,0,0,1,0,0,0,0),c(0,0,0,1,1,0,0,1)))

# may need some changes here
anova(moderator, btt=2:5)
anova(moderator, btt=7:9)
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(moderator, sigma2=1)
profile(moderator, sigma2=2)

# pub-bias indicators (will have one for each level)
level1_ability <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==3) & (moderator ==1))
lvl1abilityregtest <- regtest.rma(level1_ability)
print.regtest.rma(lvl1abilityregtest)
ranktest(level1_ability)
funnel(level1_ability)
trimfill(level1_ability, "left")
hc(level1_ability)

level2_ability <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==3) & (moderator ==2))
lvl2abilityregtest <- regtest.rma(level2_ability)
print.regtest.rma(lvl2abilityregtest)
ranktest(level2_ability)
funnel(level2_ability)
trimfill(level2_ability, "left")
hc(level2_ability)

lvl1_affect <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==1) & (moderator ==1))
lvl1affectregtest <- regtest.rma(lvl1_affect)
print.regtest.rma(lvl1affectregtest)
ranktest(lvl1_affect)
funnel(lvl1_affect)
trimfill(lvl1_affect, "right")
hc(lvl1_affect)

lvl2_affect <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==1) & (moderator ==2))
lvl2affectregtest <- regtest.rma(lvl2_affect)
print.regtest.rma(lvl2affectregtest)
ranktest(lvl2_affect)
funnel(lvl2_affect)
trimfill(lvl2_affect, "right")
hc(lvl2_affect)

lvl1_se <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==2) & (moderator ==1))
lvl1seregtest <- regtest.rma(lvl1_se)
print.regtest.rma(lvl1seregtest)
ranktest(lvl1_se)
funnel(lvl1_se)
trimfill(lvl1_se, "right")
hc(lvl1_se)

lvl2_se <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==2) & (moderator ==2))
lvl2seregtest <- regtest.rma(lvl2_se)
print.regtest.rma(lvl2seregtest)
ranktest(lvl2_se)
funnel(lvl2_se)
trimfill(lvl2_se, "right")
hc(lvl2_se)

lvl1_beh <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==4) & (moderator ==1))
lvl1behregtest <- regtest.rma(lvl1_beh)
print.regtest.rma(lvl1behregtest)
ranktest(lvl1_beh)
funnel(lvl1_beh)
trimfill(lvl1_beh, "left")
hc(lvl1_beh)

lvl2_beh <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==4) & (moderator ==2))
lvl2behregtest <- regtest.rma(lvl2_beh)
print.regtest.rma(lvl2behregtest)
ranktest(lvl2_beh)
funnel(lvl2_beh)
trimfill(lvl2_beh, "right")
hc(lvl2_beh)

lvl1_psat <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==5) & (moderator ==1))
lvl1psatregtest <- regtest.rma(lvl1_psat)
print.regtest.rma(lvl1psatregtest)
ranktest(lvl1_psat)
funnel(lvl1_psat)
trimfill(lvl1_psat)
hc(lvl1_psat)

lvl2_psat <- rma(ES_Up_Down, WT, data=reaction, subset = (DV_Code==5) & (moderator ==2))
lvl2psatregtest <- regtest.rma(lvl2_psat)
print.regtest.rma(lvl2psatregtest)
ranktest(lvl2_psat)
funnel(lvl2_psat)
trimfill(lvl2_psat, "right")
hc(lvl2_psat)

# Robust version
library(clubSandwich)
robust <- coef_test(moderator, vcov = "CR2")
robust

# Vevea & Woods
library(weightr)
abilityWoods <- reaction[reaction$DV_Code == 3,]
abknownWoods <- abilityWoods[abilityWoods$moderator==1,]
abnovelWoods <- abilityWoods[abilityWoods$moderator==2,]
abknVW2tailmod <-weightfunct(abknownWoods$ES_Up_Down, abknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abknVW2tailmod[[2]]$par
abknVW2tailsevere <-weightfunct(abknownWoods$ES_Up_Down, abknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abknVW2tailsevere[[2]]$par
abnoVW2tailmod <-weightfunct(abnovelWoods$ES_Up_Down, abnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abnoVW2tailmod[[2]]$par
abnoVW2tailsevere <-weightfunct(abnovelWoods$ES_Up_Down, abnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abnoVW2tailsevere[[2]]$par

afWoods <- reaction[reaction$DV_Code == 1,]
afknownWoods <- afWoods[afWoods$moderator==1,]
afnovelWoods <- afWoods[afWoods$moderator==2,]
afknVW2tailmod <-weightfunct(afknownWoods$ES_Up_Down, afknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
afknVW2tailmod[[2]]$par
afknVW2tailsevere <-weightfunct(afknownWoods$ES_Up_Down, afknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
afknVW2tailsevere[[2]]$par
afnoVW2tailmod <-weightfunct(afnovelWoods$ES_Up_Down, afnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
afnoVW2tailmod[[2]]$par
afnoVW2tailsevere <-weightfunct(afnovelWoods$ES_Up_Down, afnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
afnoVW2tailsevere[[2]]$par

seWoods <- reaction[reaction$DV_Code == 2,]
seknownWoods <- seWoods[seWoods$moderator==1,]
senovelWoods <- seWoods[seWoods$moderator==2,]
seknVW2tailmod <-weightfunct(seknownWoods$ES_Up_Down, seknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
seknVW2tailmod[[2]]$par
seknVW2tailsevere <-weightfunct(seknownWoods$ES_Up_Down, seknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
seknVW2tailsevere[[2]]$par
senoVW2tailmod <-weightfunct(senovelWoods$ES_Up_Down, senovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
senoVW2tailmod[[2]]$par
senoVW2tailsevere <-weightfunct(senovelWoods$ES_Up_Down, senovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
senoVW2tailsevere[[2]]$par

beWoods <- reaction[reaction$DV_Code == 4,]
beknownWoods <- beWoods[beWoods$moderator==1,]
benovelWoods <- beWoods[beWoods$moderator==2,]
beknVW2tailmod <-weightfunct(beknownWoods$ES_Up_Down, beknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
beknVW2tailmod[[2]]$par
beknVW2tailsevere <-weightfunct(beknownWoods$ES_Up_Down, beknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
beknVW2tailsevere[[2]]$par
benoVW2tailmod <-weightfunct(benovelWoods$ES_Up_Down, benovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
benoVW2tailmod[[2]]$par
benoVW2tailsevere <-weightfunct(benovelWoods$ES_Up_Down, benovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
benoVW2tailsevere[[2]]$par

psWoods <- reaction[reaction$DV_Code == 5,]
psknownWoods <- psWoods[psWoods$moderator==1,]
psnovelWoods <- psWoods[psWoods$moderator==2,]
psknVW2tailmod <-weightfunct(psknownWoods$ES_Up_Down, psknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
psknVW2tailmod[[2]]$par
psknVW2tailsevere <-weightfunct(psknownWoods$ES_Up_Down, psknownWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
psknVW2tailsevere[[2]]$par
psnoVW2tailmod <-weightfunct(psnovelWoods$ES_Up_Down, psnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
psnoVW2tailmod[[2]]$par
psnoVW2tailsevere <-weightfunct(psnovelWoods$ES_Up_Down, psnovelWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
psnoVW2tailsevere[[2]]$par

#PEESE
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==3) & (moderator==1))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==3) & (moderator==2))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==1) & (moderator==1))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==1) & (moderator==2))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==2) & (moderator==1))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==4) & (moderator==1))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==4) & (moderator==2))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==5) & (moderator==1))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (DV_Code==5) & (moderator==2))
summary(peeselm)
confint(peeselm)

