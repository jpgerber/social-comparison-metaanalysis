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
sink(file = "reaction.out", append = TRUE)
sink()
# load the dataset
reaction <- read.table("ReactionNoStapel.csv", header = TRUE, sep = ",")
covariance <- read.table("CovarianceReaction.csv", header = FALSE, sep = ",")
#for some reason covariance is read with an extra line at the end, remove it
V <- data.frame(covariance)
V <- V[-c(186),]
V<-data.matrix(V)
V <- forceSymmetric(V)
# Try it as standard random effects univariate
reaction$WT <- 1/reaction$InvWt_UpDn
overalluni <- rma(ES_Up_Down, WT, data=reaction)
hc(overalluni)
ranktest(overalluni)

# funnel, Egger's and trim-fill
funnel(overalluni)
overallunireg <- regtest.rma(overalluni)
print.regtest.rma(overallunireg)
trimfill(overalluni, "right")
# Overall (with covariance thing)
reaction$Study_Code <- relevel(factor(reaction$Study_Code), ref = "15")
reaction$Order.for.CMA <- relevel(factor(reaction$DV_Code), ref="1")

overallMV2 <- rma.mv(ES_Up_Down, V, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction)
overallMV2
# CHeck the profile plots
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(overallMV, sigma2=1)
profile(overallMV, sigma2=2)

# robust
library(robumeta)
reaction$studynum <- as.numeric(reaction$Study_Code)
overallrobu <- robu(ES_Up_Down ~ 1, data = reaction, studynum = studynum, var.eff.size = WT)

# MV with sandwich corrections
library(clubSandwich)
sandwich <- coef_test(overallMV2, vcov = "CR2")
sandwich

# then run it by PEESE
overallpeese<-rma(ES_Up_Down, WT, mods = ~I(WT), data=reaction, method = "FE")
overallpeeselm <- lm(ES_Up_Down ~ WT, data=reaction, weights = 1/WT)
summary(overallpeeselm)
confint(overallpeeselm)
# lastly, let's look at the Vevea & Woods sensitivity analysis for moderate and severe pub bias
library(weightr)
overallVW2tailmod <-weightfunct(reaction$ES_Up_Down, reaction$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .8, .6,.6,.8,.9,1))
overallVW2tailmod[[2]]$par
overallVW2tailsevere <-weightfunct(reaction$ES_Up_Down, reaction$WT, steps=c(0.05, 0.10, 0.25,0.50,.75,.9,.95,1), weights=c(1, .9, .6, .4,.4,.6,.9,1))
overallVW2tailsevere[[2]]$par