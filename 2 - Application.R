#########################################################
#            Main code to produce results in
#                                                  
#     "Objective methods for threshold assessment in 
#          weather-health watch warning systems"
#
#             Author: Pierre Masselot
#
#                    Journal
#########################################################
library(hhws)
library(splines)
library(forecast)
library(earth)
library(quantmod)
library(RColorBrewer)
library(sm)
library(devtools)
library(mgcv)
library(partykit)
library(AIM)
library(primr)
library(colorspace)

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/1 - Thresholds/Paper--Machine-Learning-Thresholds") # Repo
         
source("0 - Misc_functions.R")
source("0 - Threshold_functions.R")

#---------------------------------------------------
#                 Parameters
#---------------------------------------------------

# Parameters on data to consider
which.Y <- "Death"             # The response
which.X <- c("Tmin", "Tmax")   # The indicators
which.months <- 5:9            # the considered months

# Model parameters
Lweights <- c(.4, .4, .2)         # The lag for indicators

# Graphical parameters
cols <- c("forestgreen", "cornflowerblue", "firebrick", "goldenrod", 
  "slategrey", "black")

SAVE <- FALSE # Save the results?
result_path <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/Article V8"

#---------------------------------------------------
#                 Data Loading
#---------------------------------------------------

# Load datafile
dataread <- read.table("Data.csv", header = T, sep=";")
dataread$Date <- as.POSIXlt(apply(dataread[,1:3], 1, paste, collapse = "-"))

# Extract months of interest
datatab <- subset(dataread, Month %in% which.months)

# Extract variables of interest
datatab <- datatab[,c("Date", which.Y, which.X)]

# Number of years, and days of data, and number of weather variables
nyear <- length(unique(datatab$Date$year))
n <- nrow(datatab)
p <- length(which.X)

#---------------------------------------------------
#             Compute Over-mortality
#---------------------------------------------------

# Compute baseline through natural splines
em.form <- sprintf("%s ~ ns(Date$yday, 4) + ns(Date$year, round(nyear / 10))",
  which.Y)
EM <- lm(as.formula(em.form), datatab)$fitted.values

# Compute over-mortality
datatab$OM <- excess(datatab[,which.Y], EM)

#---------------------------------------------------
#               Prepare indicators
#---------------------------------------------------

L <- length(Lweights)
# indicator.names <- sprintf("Z%s", which.X)
indicator.names <- which.X

# Compute the moving-average
indicators <- apply(datatab[,which.X], 2, filter, 
  Lweights, sides = 1)
# To account for breaks in the series
indicators[which(diff(datatab$Date, L) > L) + L,] <- NA 

# Adding indicators to the data table
datatab[,indicator.names] <- indicators

#---------------------------------------------------
#            Remove NAs from data table
#---------------------------------------------------

datatab <- na.omit(datatab)
n <- nrow(datatab)

#---------------------------------------------------
#            Algorithms application
#---------------------------------------------------
results <- list()

#----- Model-based regression trees (MOB) -----
treefit <- lmtree(OM ~ Tmin + Tmax, data = datatab, minsize = 4)

# Extract thresholds from the tree object
results[["MOB"]] <- extractThresholds_party(treefit)$thresholds

# Plot the tree
x11(width = 17, height = 7)
plot(treefit, ip_args = list(id = FALSE, pval = FALSE), 
  tp_args = list(id = FALSE))
  
if (SAVE){
  dev.print(png, filename = sprintf("%s/FigS1_MOB.png", result_path), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/FigS1_MOB.eps", result_path))
}

#----- Multivariate adaptive regression splines (MARS) -----
marsfit <- earth(OM ~ Tmin + Tmax, data = datatab, degree = p, endspan = 4)

# Extract thresholds
results[["MARS"]] <- apply(marsfit$cuts[,indicator.names], 2, max)

# Plot the surface
surfgrid <- do.call(expand.grid, 
  lapply(datatab[,indicator.names], 
    function(x) seq(min(x), max(x), length.out = 1000) )
)
surf <- matrix(predict(marsfit, surfgrid, type = "response"), nrow = 1000)

x11()
filled.contour(unique(surfgrid[,1]), unique(surfgrid[,2]), surf,
  color.palette = diverge_hcl, zlim = c(-110, 110), cex.axis = 1.2,
  key.title = title(main = "OM", cex = 1.3), 
  plot.title = {
    title(xlab = indicator.names[1], cex.lab = 1.3)
    title(ylab = indicator.names[2], cex.lab = 1.3)
    abline(v = unique(marsfit$cuts[marsfit$cuts[,1] != 0,1]), lty = 2)
    abline(v = results[["MARS"]][1], lwd = 2)
    mtext(formatC(results[["MARS"]][1], format = "f", digits = 1), 
      at = results[["MARS"]][1], line = 0.5, cex = 1.5)
    abline(h = unique(marsfit$cuts[marsfit$cuts[,2] != 0,2]), lty = 2)
    abline(h = results[["MARS"]][2], lwd = 2)
    mtext(formatC(results[["MARS"]][2], format = "f", digits = 1), 
      at = results[["MARS"]][2], line = 0.5, cex = 1.5, side = 2)
  }
)

if (SAVE){
  dev.print(png, filename = sprintf("%s/FigS2_MARS.png", result_path), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/FigS2_MARS.eps", result_path))
}


#----- Patient rule-induction method (PRIM) -----
# Initial peeling
peelres <- peeling(datatab$OM, datatab[,indicator.names], 
  beta.stop = 4/n, peeling.side = -1)

# Plot and analyze the trajectory
x11()
plot_trajectory(peelres, pch = 16, col = cols[3], ylab = "Mean OM", 
  xlim = c(0, 0.02), support = 6/n, abline.pars = list(lwd = 2, lty = 2))
mtext("n = 6", at = 6/n, cex = 1.5, line = 0.5)

if (SAVE){
  dev.print(png, filename = sprintf("%s/FigS3_PRIM.png", result_path), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/FigS3_PRIM.eps", result_path))
}

# Thresholds refinement by pasting 
primres <- pasting(peelres, support = 6/n)

# Thresholds extraction
results[["PRIM"]] <- sapply(extract.box(primres)$limits[[1]], "[", 1) 

#----- Adaptive index models (AIM) -----
# Cross-validation to obtain the optimal number of steps
cvb <- cv.lm.main(datatab[,indicator.names], datatab$OM, 
  nsteps = 3 * p, maxnumcut = 3, backfit = T, mincut = 4/n)

# Fit the AIM with the optimal number of steps
aimfit <- lm.main(datatab[,indicator.names], datatab$OM, nsteps = cvb$kmax,
  maxnumcut = 3, backfit = T, mincut = 4/n)

# Extract thresholds
rules <- aimfit$res[[cvb$kmax]]
maxrules <- aggregate(rules[,2], by = list(var = rules[,1]), max)
results[["AIM"]] <- rep(NA_real_, p)
results[["AIM"]][maxrules$var] <- maxrules$x

# Plot the index
surfaim <- matrix(index.prediction(rules, surfgrid), nrow = 1000)

x11()
filled.contour(unique(surfgrid[,1]), unique(surfgrid[,2]), surfaim,
  color.palette = function(n) hcl.colors(n, "YlOrBr", rev = T), cex.axis = 1.2,
  key.title = title(main = "Score", cex = 1.3), levels = -1:cvb$kmax + 0.5,
  key.axes = axis(4, at = 0:cvb$kmax + 0.5, labels = 0:cvb$kmax),
  plot.title = {
    title(xlab = indicator.names[1], cex.lab = 1.3)
    title(ylab = indicator.names[2], cex.lab = 1.3)
    if (!is.na(results[["AIM"]][1])){
      abline(v = rules[rules[,1] == 1,2], lty = 2)
      abline(v = maxrules$x[maxrules$var == 1], lwd = 2)
      mtext(formatC(results[["AIM"]][1], format = "f", digits = 1), 
        at = results[["AIM"]][1], line = 0.5, cex = 1.5)
    }
    if (!is.na(results[["AIM"]][2])){
      abline(h = rules[rules[,1] == 2,2], lty = 2)
      abline(h = maxrules$x[maxrules$var == 2], lwd = 2)
      mtext(formatC(results[["AIM"]][2], format = "f", digits = 1), 
        at = results[["AIM"]][2], line = 0.5, cex = 1.5, side = 2)
    }
  }
)

if (SAVE){
  dev.print(png, filename = sprintf("%s/FigS4_AIM.png", result_path), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/FigS4_AIM.eps", result_path))
}


#----- Generalized additive models (GAM) -----
gamfit <- gam(OM ~ s(Tmin) + s(Tmax), data = datatab)

# Extract thresholds
results[["GAM"]] <- extractThresholds_gam(gamfit)

# Plot the fitted functions with thresholds
x11(width = 15)
par(mfrow = c(1,2))
for (j in 1:p){
  plot(gamfit, select = j, shade = TRUE, shade.col = cols[5], lwd = 2, 
    rug = FALSE, cex.lab = 1.3, cex.axis = 1.2, ylab = "Estimated function")
  abline(h = 0, lty = 3)
  abline(v = results[["GAM"]][j], lwd = 2, lty = 2)
  mtext(formatC(results[["GAM"]][j]), at = results[["GAM"]][j], 
    cex = 1.5, line = 0.5)
}

if (SAVE){
  dev.print(png, filename = sprintf("%s/FigS5_GAM.png", result_path), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/FigS5_GAM.eps", result_path))
}

# Saving estimated thresholds
if (SAVE){
  thresholds <- Reduce(rbind, results)
  rownames(thresholds) <- names(results)
  write.table(thresholds, 
    file = sprintf("%s/Results_application.csv", result_path),
    quote = FALSE, sep = ";")
}

#---------------------------------------------------
#            Results analysis
#---------------------------------------------------

# Add reference (Chebana et al. 2013)
results[["Reference"]] <- c(20, 33)

# Extract alarms for each method
alarms <- lapply(results, extract_alarms, x = datatab[,indicator.names], 
  y = datatab$OM)

# Analyse basic statistics of alarms
alarms.n <- sapply(alarms, length)
episodes.n <- sapply(alarms, 
  function(x) sum(diff(as.integer(names(x))) > 3) + 1)
alarms.mean <- sapply(alarms, mean)
alarms.coverage <- alarms.mean * alarms.n / n

# Compare to a predefined cutpoint
cutPts <- seq(30, 60, 5)
sensitivity <- falseAlarms <- specificity <- precision <- F1 <- F2 <- list()
for (i in 1:length(cutPts)){
  trueAlarms <- episodes(datatab$OM, cutPts[i])
  found <- lapply(alarms, function(x) trueAlarms[trueAlarms$t %in% names(x),])
  sensitivity[[i]] <- sapply(found, function(x) nrow(x) / nrow(trueAlarms))
  falseAlarms[[i]] <- sapply(alarms, 
    function(x) sum(!names(x) %in% trueAlarms$t))
  specificity[[i]] <- 1 - falseAlarms[[i]] / (n - nrow(trueAlarms))
  precision[[i]] <- sapply(alarms, function(x) 
    mean(as.integer(names(x)) %in% trueAlarms$t))
  F1[[i]] <- sapply(alarms, function(x) 
    Fscore(as.integer(names(x)), trueAlarms$t, 1))
  F2[[i]] <- sapply(alarms, function(x) 
    Fscore(as.integer(names(x)), trueAlarms$t, 2))
}

#---- Plot results ----
x11(height = 10)
par(mfrow = c(4, 1), mar = c(5, 5, 4, 2) + .1)
matplot(cutPts, t(as.data.frame(sensitivity)), type = "b", col = cols,
  pch = 14 + 1:length(results), ylab = "Sensitivity", 
  cex = 1.5, lwd = 3, cex.lab = 1.5, cex.axis = 1.2, xlab = "")
text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "a)", 
  cex = 3, xpd = T, adj = c(-0.5,1.2))
outerLegend("topcenter", names(results), lwd = 3, pch = 14 + 1:length(results), 
  lty = 1:(length(results) + 1), col = cols, ncol = 4, bty = "n", cex = 1.5,
  pt.cex = 2)
matplot(cutPts, t(as.data.frame(precision)), type = "b", col = cols,
  pch = 14 + 1:length(results),  ylab = "Precision", 
  cex = 1.5, lwd = 3, cex.lab = 1.5, cex.axis = 1.2, xlab = "")
text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "b)", 
  cex = 3, xpd = T, adj = c(-0.5,1.2))
matplot(cutPts, t(as.data.frame(F1)), type = "b", col = cols,
  pch = 14 + 1:length(results),  ylab = expression(F[1]), 
  cex = 1.5, lwd = 3, cex.lab = 1.5, cex.axis = 1.2, xlab = "")
text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "c)", 
  cex = 3, xpd = T, adj = c(-0.5,1.2))
matplot(cutPts, t(as.data.frame(F2)), type = "b", col = cols,
  pch = 14 + 1:length(results), xlab = "Cut point", ylab = expression(F[2]), 
  cex = 1.5, lwd = 3, cex.lab = 1.5, cex.axis = 1.2)
text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "d)", 
  cex = 3, xpd = T, adj = c(-0.5,1.2))
  
if (SAVE){
  dev.print(png, filename = sprintf("%s/Fig4_Performance.png", result_path), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/Fig4_Performance.eps", result_path))
}



#---------------------------------------------------
#            Bootstrap simulations
#---------------------------------------------------
B <- 1000

# Create data blocks
datablock <- split(datatab, datatab$Date$year)
bpool <- 1:length(datablock)

# Draw bootstrap samples
bsamples <- replicate(B, sample(datablock, replace = T), simplify = F)
bsamples <- lapply(bsamples, Reduce, f = "rbind")

# Apply algorithms
bootRes <- list()
bootRes[["MOB"]] <- lapply(bsamples, function(dat) 
  MOB.apply(dat$OM, dat[,indicator.names],  minsize = 4))
bootRes[["MARS"]] <- lapply(bsamples, function(dat) 
  MARS.apply(dat$OM, dat[,indicator.names], degree = p, endspan = 4))
bootRes[["PRIM"]] <- lapply(bsamples, function(dat){ 
  PRIM.apply(dat$OM, dat[,indicator.names], beta.stop = 4/n, peeling.side = -1)  
})
bootRes[["AIM"]] <- lapply(bsamples, function(dat) 
  AIM.apply(dat$OM, dat[,indicator.names], backfit = T, mincut = 4/n))
bootRes[["GAM"]] <- lapply(bsamples, function(dat) 
  GAM.apply(dat$OM, dat[,indicator.names]))
 
if (SAVE){
  save(bsamples, bootRes, 
    file = sprintf("%s/Results_bootstrap.RData", result_path))
}

# Obtain thresholds
bootThresh <- lapply(bootRes, sapply, "[[", "thresholds")

#---- Plot the bootstrap results
nm <- length(results) - 1
ylim <- range(sapply(bootRes, sapply, "[[", "thresholds"), na.rm = TRUE)
xlim <- c(.5, p * nm + 0.5)

x11(width = 10, height = 6)
par(yaxs = "i", cex.axis = 1.2, cex.lab = 1.3, mar = c(5, 2, 4, 1) + .1)
layout(matrix(1:2, nrow = 1), width = c(0.7, 0.3))
for (j in 1:p){    
  resj <- sapply(bootThresh, "[", j, )
  boxplot(resj, horizontal = T, at = (1:nm) * p - j + 1, yaxt = "n",
    border = cols, varwidth = T, 
    lwd = 2, ylim = ylim, xlim = xlim, add = j > 1, pch = j, 
    xlab = "Temperature", main = "a) Bootstrap thresholds")
  obsThresh <- sapply(results[1:nm], "[", j)
  segments(obsThresh, (1:nm - 1) * p + 0.48, 
    obsThresh, (1:nm) * p + 0.48, 
    lwd = 3, lty = 1)
}
abline(h = 0.5 + p * 1:(nm - 1), lty = 3)
axis(side = 2, at = (1:nm) * p - (p - 1)/2, labels = names(bootRes),
  tick = F)

# Number of unfound thresholds
percentNA <- sapply(bootThresh, apply, 1, function(x) mean(is.na(x)))
bp <- barplot(percentNA[p:1,], beside = TRUE, horiz = TRUE, xlim = c(0,1),
  col = rep(cols[1:nm], each = p), space = c(0, 0.1), xaxt = "n",
  xlab = "Proportion (%)", main = "b) Failure proportion")
axis(1, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2) * 100)
abline(h = bp[1, -1] - 0.55, lty = 3)
box()

if (SAVE){
  dev.print(png, filename = sprintf("%s/Fig3_Bootstrap.png", result_path), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/Fig3_Bootstrap.eps", result_path))
}


