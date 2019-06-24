
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
reaction$moderator <- reaction$SimDissim
# Univariate approach
reaction$moderator <- relevel(factor(reaction$moderator), ref="1")
reaction$DV_Code <- relevel(factor(reaction$DV_Code),ref="3")
uni <- rma(ES_Up_Down, WT, mods = ~ DV_Code*moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
uni
anova(uni, btt=2:5)
anova(uni, btt=7:8)
reviseduni <- rma(ES_Up_Down, WT, mods = ~ moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
reviseduni
addmargins(table(reaction$moderator))
predict(reviseduni, newmods = rbind(c(0),c(1)))


# Multivariate approach
reaction$moderator <- relevel(factor(reaction$moderator), ref="1")
reaction$DV_Code <- relevel(factor(reaction$DV_Code),ref="3")
moderator <- rma.mv(ES_Up_Down, V, mods = ~ DV_Code*moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
addmargins(table(reaction$DV_Code, reaction$moderator))
moderator
# may need some changes here
anova(moderator, btt=2:5)
anova(moderator, btt=7:8)
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(moderator, sigma2=1)
profile(moderator, sigma2=2)

# Try a reduced model (no interaction or DV code)
moderatorreduced <- rma.mv(ES_Up_Down, V, mods = ~ moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
moderatorreduced
profile(moderatorreduced, sigma2=1)
profile(moderatorreduced, sigma2=2)
moderatorreducedmore <- rma.mv(ES_Up_Down, V, mods = ~ moderator, random = ~ 1 | Order.for.CMA, data=reaction, knha=TRUE)
moderatorreducedmore
profile(moderatorreducedmore, sigma2=1)
predict(moderatorreducedmore, newmods = rbind(c(0),c(1)))

# pub-bias indicators (will have one for each level)
simtrim <- rma(ES_Up_Down, WT, data=reaction, subset = (moderator ==1))
simtrim
simtrimregtest <- regtest.rma(simtrim)
print.regtest.rma(simtrimregtest)
trimfill(simtrim, "left")
funnel(simtrim)
ranktest(simtrim)
hc(simtrim)

dissimtrim <- rma(ES_Up_Down, WT, data=reaction, subset = (moderator ==2))
dissimtrimregtest <- regtest.rma(dissimtrim)
print.regtest.rma(dissimtrimregtest)
trimfill(dissimtrim, "right")
funnel(dissimtrim)
ranktest(dissimtrim)
hc(dissimtrim)
# Robust version
library(clubSandwich)
simdissim_robust <- coef_test(moderatorreducedmore, vcov = "CR2")
simdissim_robust

# Vevea & Woods
library(weightr)
simWoods <- reaction[reaction$moderator==2,]
dissimWoods <- reaction[reaction$moderator==1,]
abknVW2tailmod <-weightfunct(simWoods$ES_Up_Down, simWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abknVW2tailmod[[2]]$par
abknVW2tailsevere <-weightfunct(simWoods$ES_Up_Down, simWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abknVW2tailsevere[[2]]$par
abnoVW2tailmod <-weightfunct(dissimWoods$ES_Up_Down, dissimWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abnoVW2tailmod[[2]]$par
abnoVW2tailsevere <-weightfunct(dissimWoods$ES_Up_Down, dissimWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abnoVW2tailsevere[[2]]$par

#PEESE

peese1lm <- lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (moderator==1))
summary(peese1lm)
confint(peese1lm)
peese2lm <- lm(ES_Up_Down ~WT, data =reaction,weight=1/WT, subset = (moderator==2))
summary(peese2lm)
confint(peese2lm)
