# clear workspace
rm(list=ls())
# store the current directory
initial.dir<-getwd()
# change to the new directory
setwd("/Users/jonathan.gerber/Dropbox/SCMeta_Articles/000 SCMeta Reports/7 Second resubmit/Analyses")
# load the necessary libraries
library(metafor)
# load the dataset
selection <- read.table("SCMeta_SelectionNOsimREAL.csv", header = TRUE, sep = ",")
library(lme4)
# overall effect
attach(selection)
selection$xi <- round(xi,0)
selection$ni <- round(ni,0)
selection$pi <- selection$xi / selection$ni
selection$plox <- log(selection$pi / (1-selection$pi))
detach(selection)

table1 <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=selection)
table1
funnel(table1, yaxis = "ni")
regtesttable1 <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ ni, data=selection)
regtesttable1
regtesttable2 <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ ni, data=selection, subset = (dimension==1) & (setting==1))
regtesttable2
regtesttable3 <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ ni, data=selection, subset = ((dimension==3) & (setting==1)))
regtesttable3
regtesttable4 <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ ni, data=selection, subset = (dimension==1) & (setting==5))
regtesttable4
regtesttable5 <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ ni, data=selection, subset = (dimension==2) & (setting==5))
regtesttable5
regtesttable6 <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ ni, data=selection, subset = (dimension==1) & (setting==6))
regtesttable6

# Now for pub bias do it as the other version
selectionESs <- escalc(measure="PLO", xi=xi, ni=ni, data=selection)
res <- rma(yi, vi, data=selectionESs)
selectionESs$ninew <- selection$ni
res <- rma(yi, 1/ninew, data=selectionESs)

### check for funnel plot asymmetry
funnel(res)
regtest(res)

### regression test is just the same as using sqrt(vi) as moderator
res <- rma(yi, vi, mods = ~ sqrt(vi), data=dat)
res
#Now with both moderators
selection$dimension <- relevel(factor(selection$dimension), ref="1")
selection$setting <- relevel(factor(selection$setting), ref="1")
addmargins(table(selection$dimension, selection$setting))

doublemod <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ setting*dimension,  data=selection)
doublemod
predict(doublemod, newmods=c(0,0,0,0))
predict(doublemod, newmods=c(1,0,0,0))
predict(doublemod, newmods=c(0,1,0,0))
predict(doublemod, newmods=c(1,0,1,0))
predict(doublemod, newmods=c(0,0,0,1))

predict(doublemod, transf=transf.ilogit)

# then we need to go over the near/far stuff
# Getting the right subtotals / variables
attach(selection)
selection$nearxi <- round((nearperc/100 *ni),0)
selection$farxi <- round((farperc/100 *ni),0)
selection$nfni <- farxi + nearxi
selection$nearbixi <- round((nearbi/100 *ni),0)
selection$farbixi <- round((farbi/100 *ni),0)
selection$nfbini <- farbixi + nearbixi
detach(selection)

###### one direction
rma.glmm(measure="PLO", xi=nearxi, ni=nfni, data=selection)

# Then with threat (up version)
threatselectionnearfar <- rma.glmm(measure="PLO", xi=nearxi, ni=nfni, mods = ~ dimension,  data=selection)
threatselectionnearfar

selectionup <- selection[selection$nfni > 0,]
addmargins(table(selectionup$dimension))

# try the two mod
doublemodnearfar <- rma.glmm(measure="PLO", xi=nearxi, ni=nfni, mods = ~ setting*dimension,  data=selection)
doublemodnearfar
# Then calculate the downward variable.
attach(selection)
selection$neardownxi <- nearbixi - nearxi
selection$fardownxi <- farbixi - farxi
selection$downni <- neardownxi + fardownxi
detach(selection)
threatselectiondownnearfar <- rma.glmm(measure="PLO", xi=neardownxi, ni=downni, mods = ~ dimension,  data=selection)
threatselectiondownnearfar

#run bi-directional
rma.glmm(measure="PLO", xi=nearbixi, ni=nfbini, data=selection)
# with threat
threatselectionnearfarbi <- rma.glmm(measure="PLO", xi=nearbixi, ni=nfbini, mods = ~ dimension,  data=selection)
threatselectionnearfarbi
# try the two mod
doublemodnearfarbi <- rma.glmm(measure="PLO", xi=nearbixi, ni=nfbini, mods = ~ setting*dimension,  data=selection)
doublemodnearfarbi
# close the output file
sink()

### NOw here's the tri-chotomous stuff
# First, create the file with three different subtotals for comparison and the correct xis.
selectionsim <- read.table("SCMeta_SelectionsimLONG.csv", header = TRUE, sep = ",")

attach(selectionsim)
selectionsim$xi <- round(xi,0)
selectionsim$ni <- round(ni,0)
detach(selectionsim)

selectionsim$dimension <- relevel(factor(selectionsim$dimension), ref="1")
selectionsim$setting <- relevel(factor(selectionsim$setting), ref="1")
addmargins(table(selectionsim$dimension, selectionsim$setting))

# Up vs sim (xi vs zi)
upneutraldata <- read.table("triupsim.csv", header = TRUE, sep = ",")
updowndata <- read.table("triupdown.csv", header = TRUE, sep = ",")
downneutraldata <- read.table("tridownsim.csv", header = TRUE, sep = ",")

upneutral <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods =~ direction, data=upneutraldata)
upneutral
updown <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods =~ direction, data=updowndata)
updown
downsim <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods =~ direction, data=downneutraldata)
downsim
upneutralmods <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods =~ direction*dimension*setting, data=upneutraldata)
upneutralmods
updownmods <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods =~ dimension*direction*setting, data=updowndata)
updownmods
downsimmods <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods =~ dimension*direction*setting, data=downneutraldata)
downsimmods
# Make the averages

up_av <- sum(selectionsim$xi)/sum(selectionsim$ni)
up_av
down_av <- sum(selectionsim$yi)/sum(selectionsim$ni)
down_av
sim_av <- sum(selectionsim$zi)/sum(selectionsim$ni)
sim_av

down_av <- mean(selectionsim$downperc)
sim_av <- mean(selectionsim$simperc)


############ Now trying some pub bias stuff ############### 

# unload the libraries
detach("package:metage")
# change back to the original directory
setwd(initial.dir)