
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

############# Objective Subjective DV ################
# Set the moderator name
reaction$moderator <- reaction$Obj_Subj_REAL
# Univariate approach
reaction$moderator <- relevel(factor(reaction$moderator), ref="1")
reaction$DV_Code <- relevel(factor(reaction$DV_Code),ref="3")
uni <- rma(ES_Up_Down, WT, mods = ~ DV_Code*moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
uni
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
# may need some changes here
anova(moderator, btt=2:5)
anova(moderator, btt=7:9)
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(moderator, sigma2=1)
profile(moderator, sigma2=2)

######### SELF-OTHER PART ##############
#Set the moderator name
reaction$moderator <- reaction$SelfOtherDirection
# Univariate approach
reaction$moderator <- relevel(factor(reaction$moderator), ref="1")
reaction$DV_Code <- relevel(factor(reaction$DV_Code),ref="3")
uni <- rma(ES_Up_Down, WT, mods = ~ DV_Code*moderator, random = ~ 1 | Study_Code/Order.for.CMA, data=reaction, knha=TRUE)
uni
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
# may need some changes here
anova(moderator, btt=2:5)
anova(moderator, btt=7:9)
par(mar=c(1,1,1,1))
par(mfrow=c(2,1))
profile(moderator, sigma2=1)
profile(moderator, sigma2=2)
