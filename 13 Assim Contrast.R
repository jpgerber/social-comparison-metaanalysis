# Draft of control groups section.
rm(list=ls())
initial.dir<-getwd()
setwd("/Users/jonathan.gerber/Dropbox/SCMeta_Articles/000 SCMeta Reports/7 Second resubmit/Analyses")
library(metafor)
assimcontrast <- read.table("AssimContrastLong.csv", header = TRUE, sep = ",")
assimcontrastcovariance <- read.table("AssimContrastCov.csv", header = FALSE, sep = ",")

assimcontrast$WT <- 1/assimcontrast$InvWt

attach(assimcontrast)
assimcontrast$Study_Code <- round(Study_Code,0)
detach(assimcontrast)

V <- data.frame(assimcontrastcovariance)
V <- V[-c(186),]
V<-data.matrix(V)
V <- forceSymmetric(V, uplo = "L")

# Overall check of whether DV does something (univariate)
assimcontrast$DV_Code <- relevel(factor(assimcontrast$DV_Code),ref="3")
byDV <- rma.mv(ES, WT, mods = ~ DV_Code, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
byDV
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(byDV, sigma2=1)
profile(byDV, sigma2=2)

# Multivariate check of DV
assimcontrast$DV_Code <- relevel(factor(assimcontrast$DV_Code),ref="3")
byDV <- rma.mv(ES, V, mods = ~ DV_Code, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
byDV
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(byDV, sigma2=1)
profile(byDV, sigma2=2)

# Looking at just up and down (univariate)
assimcontrast$Index1 <- relevel(factor(assimcontrast$Index1),ref="1")
Direction <- rma.mv(ES, WT, mods = ~ Index1, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
Direction
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(Direction, sigma2=1)
profile(Direction, sigma2=2)
predict(Direction, newmods = rbind(c(0),c(1)))

# Multivariate directional effects
assimcontrast$Index1 <- relevel(factor(assimcontrast$Index1),ref="1")
Direction <- rma.mv(ES, V, mods = ~ Index1, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
Direction
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(Direction, sigma2=1)
profile(Direction, sigma2=2)
predict(Direction, newmods = rbind(c(0),c(1)))

# Pub bias for the directional effect
Directionuni <- rma(ES, WT, mods = ~ Index1, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
Directionuni

up <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==1))
upregtest <- regtest.rma(up)
print.regtest.rma(upregtest)
ranktest(up)
funnel(up)
trimfill(up, "left")
hc(up)

down <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==2))
downregtest <- regtest.rma(down)
print.regtest.rma(downregtest)
ranktest(down)
funnel(down)
trimfill(down, "left")
hc(down)
#RVE
library(clubSandwich)
robust <- coef_test(Direction, vcov = "CR2")
robust

# Vevea & Woods
library(weightr)
upWoods <- assimcontrast[assimcontrast$Index1==1,]
downWoods <- assimcontrast[assimcontrast$Index1==2,]
abknVW2tailmod <-weightfunct(upWoods$ES, upWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abknVW2tailmod[[2]]$par
abknVW2tailsevere <-weightfunct(upWoods$ES, upWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abknVW2tailsevere[[2]]$par
abnoVW2tailmod <-weightfunct(downWoods$ES, downWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
abnoVW2tailmod[[2]]$par
abnoVW2tailsevere <-weightfunct(downWoods$ES, downWoods$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
abnoVW2tailsevere[[2]]$par

#Peese
uppeese<-rma(ES, WT, mods = ~I(WT), data=assimcontrast, method = "FE", subset = (Index1==1))
uppeese
downpeese<-rma(ES, WT, mods = ~I(WT), data=assimcontrast, method = "FE", subset = (Index1==2))
downpeese
up.peese=lm(ES ~WT, data =upWoods,weight=1/WT)
summary(up.peese)
confint(up.peese)
down.peese=lm(ES ~WT, data =downWoods,weight=1/WT)
summary(down.peese)
confint(down.peese)

# Then run a very large meta-regression model
# Have to drop SimDissim, SelfOther and Obj_Subj as too few values. (univariate)
LargeModel <- rma.mv(ES, WT, mods = ~ Task*Index1 + TargetDistance*Index1 
                     + Knowndimension*Index1 + Whovaried*Index1 + Exemplar*Index1, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
LargeModel
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(LargeModel, sigma2=1)
profile(LargeModel, sigma2=2)

# Multivariate large model.
MVLargeModel <- rma.mv(ES, V, mods = ~ Task*Index1 + TargetDistance*Index1 
                     + Knowndimension*Index1 + Whovaried*Index1 + Exemplar*Index1, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
MVLargeModel
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(MVLargeModel, sigma2=1)
profile(MVLargeModel, sigma2=2)

# Then drop the non-sig predictors (multivariate)
ReducedModel1 <- rma.mv(ES, V, mods = ~ TargetDistance*Index1 
                     + Whovaried*Index1, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
ReducedModel1
profile(ReducedModel1, sigma2=1)
profile(ReducedModel1, sigma2=2)
# Getting cell means & CIs
predict(ReducedModel1, newmods = rbind(c(0,0,0,0,0),c(0,1,0,0,0),c(1,0,0,0,0),c(1,1,0,1,0),c(0,0,1,0,0),c(0,1,1,0,1),
                                       c(1,0,1,0,0),c(1,1,1,1,1)))

ReducedModel2 <- rma.mv(ES, V, mods = ~ TargetDistance*Index1 
                        + Whovaried*Index1, random = ~ 1 | id, data=assimcontrast, knha=TRUE)
ReducedModel2
profile(ReducedModel1, sigma2=1)
# Getting cell means & CIs
predict(ReducedModel1, newmods = rbind(c(0,0,0,0,0),c(0,1,0,0,0),c(1,0,0,0,0),c(1,1,0,1,0),c(0,0,1,0,0),c(0,1,1,0,1),
                                       c(1,0,1,0,0),c(1,1,1,1,1)))
addmargins(table(assimcontrast$TargetDistance, assimcontrast$Index1, assimcontrast$Whovaried))
# MV reduced model (THIS IS TOO SMALL!)
MVReducedModel1 <- rma.mv(ES, V, mods = ~ Whovaried*Index1, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
MVReducedModel1
profile(MVReducedModel1, sigma2=1)
profile(MVReducedModel1, sigma2=2)
predict(MVReducedModel1, newmods = rbind(c(0,0,0),c(0,1,0),c(1,0,0),c(1,1,1)))

# Also run the highly reduced model of sim/dissim (univariate)
SimDissimModel <- rma.mv(ES, WT, mods = ~ Index1*SimDissim, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
SimDissimModel
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(SimDissimModel, sigma2=1)
profile(SimDissimModel, sigma2=2)

# Sim/dissim (MULTIVARIATE)
MVSimDissimModel <- rma.mv(ES, V, mods = ~ Index1*SimDissim, random = ~ 1 | Study_Code/id, data=assimcontrast, knha=TRUE)
MVSimDissimModel
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(MVSimDissimModel, sigma2=1)
profile(MVSimDissimModel, sigma2=2)

#### And here's the pub bias stuff for the three-way interaction
uplocalself <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==1) & (TargetDistance ==1) & (Whovaried ==1))
trimregtest <- regtest.rma(uplocalself)
print.regtest.rma(trimregtest)
trimfill(uplocalself, "left")
funnel(uplocalself)
ranktest(uplocalself)

downlocalself <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==2) & (TargetDistance ==1) & (Whovaried ==1))
trimregtest <- regtest.rma(downlocalself)
print.regtest.rma(trimregtest)
trimfill(downlocalself, "left")
funnel(downlocalself)
ranktest(downlocalself)

updistself <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==1) & (TargetDistance ==1) & (Whovaried ==2))
trimregtest <- regtest.rma(uplocalstand)
print.regtest.rma(trimregtest)
trimfill(uplocalstand, "left")
funnel(uplocalstand)
ranktest(uplocalstand)

downdistself <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==2) & (TargetDistance ==1) & (Whovaried ==2))
trimregtest <- regtest.rma(downdistself)
print.regtest.rma(trimregtest)
trimfill(downdistself, "left")
funnel(downdistself)
ranktest(downdistself)

uplocalstand <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==1) & (TargetDistance ==2) & (Whovaried ==1))
trimregtest <- regtest.rma(uplocalstand)
print.regtest.rma(trimregtest)
trimfill(uplocalstand, "left")
funnel(uplocalstand)
ranktest(uplocalstand)

downlocalstd <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==2) & (TargetDistance ==2) & (Whovaried ==1))
trimregtest <- regtest.rma(downlocalstd)
print.regtest.rma(trimregtest)
trimfill(downlocalstd, "left")
funnel(downlocalstd)
ranktest(downlocalstd)

updiststand <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==1) & (TargetDistance ==2) & (Whovaried ==2))
trimregtest <- regtest.rma(updiststand)
print.regtest.rma(trimregtest)
trimfill(updiststand, "left")
funnel(updiststand)
ranktest(updiststand)

downdiststand <- rma(ES, WT, data=assimcontrast, subset = (Index1 ==2) & (TargetDistance ==2) & (Whovaried ==2))
trimregtest <- regtest.rma(downdiststand)
print.regtest.rma(trimregtest)
trimfill(downdiststand, "left")
funnel(downdiststand)
ranktest(downdiststand)

library(weightr)
up <- assimcontrast[assimcontrast$Index1 == 1,]
down <- assimcontrast[assimcontrast$Index1 == 2,]
uplocal <- up[up$TargetDistance == 1,]
updist <- up[up$TargetDistance == 2,]
downlocal <- down[down$TargetDistance == 1,]
downdist <- down[down$TargetDistance == 2,]
eight1way <- uplocal[uplocal$Whovaried ==1,]
eight2way <- uplocal[uplocal$Whovaried ==2,]
eight3way <- updist[updist$Whovaried ==1,]
eight4way <- updist[updist$Whovaried ==2,]
eight5way <- downlocal[downlocal$Whovaried ==1,]
eight6way <- downlocal[downlocal$Whovaried ==2,]
eight7way <- downdist[downdist$Whovaried ==1,]
eight8way <- downdist[downdist$Whovaried ==2,]

# I just change it here to save the lines.
Woodsdata <- eight8way
moderate <-weightfunct(Woodsdata$ES, Woodsdata$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
moderate[[2]]$par
severe <-weightfunct(Woodsdata$ES, Woodsdata$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
severe[[2]]$par

#PEESE
peeselm <-lm(ES ~WT, data =assimcontrast,weight=1/WT, subset = (Index1 ==1) & (Whovaried==1) & (TargetDistance ==1))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES ~WT, data =assimcontrast,weight=1/WT, subset = (Index1 ==2) & (Whovaried==1) & (TargetDistance ==1))
summary(peeselm)
confint(peeselm)
peeselm <-lm(ES ~WT, data =assimcontrast,weight=1/WT)
summary(peeselm)
confint(peeselm)
