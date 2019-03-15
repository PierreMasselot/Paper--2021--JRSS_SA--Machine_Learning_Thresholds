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
library(rpart)
library(splines)
library(forecast)
library(earth)
library(quantmod)
library(RColorBrewer)
library(sm)

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/1 - Thresholds/Paper--Machine-Learning-Thresholds") # Data path
              
source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/PRIM/R/PRIM_functions.R")
source("Functions_rpartExtractResults.R")
source("Threshold_functions.R")
source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Plot_functions.R")
source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Marc/Plot_functions_Marc.R")
source("Other_functions.R")

#---------------------------------------------------
#                 Parameters
#---------------------------------------------------

# Parameters on data to consider
which.Y <- "Death"             # The response
which.X <- c("Tmin", "Tmax")   # The indicators
which.months <- 5:9            # the considered months

# Model parameters
Lweights <- c(.4, .4, .2)         # The lag for indicators

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
indicator.names <- sprintf("Z%s", which.X)

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

#---------------------------------------------------
#            Algorithms application
#---------------------------------------------------

results <- list()
results[["CART"]] <- CART.apply(datatab$OM, datatab[,indicator.names])
results[["MARS"]] <- MARS.apply(datatab$OM, datatab[,indicator.names], p = p)
results[["PRIM"]] <- PRIM.apply(datatab$OM, datatab[,indicator.names])

#---------------------------------------------------
#            Results analysis
#---------------------------------------------------

# Add reference (Chebana et al. 2013)
refThresh <- c(20, 33)
uni.alb <- mapply(">=", datatab[,indicator.names], refThresh)
alarmlogi <- apply(uni.alb, 1, all)
alarmRef <- datatab$OM[alarmlogi]
names(alarmRef) <- which(alarmlogi)
results[["Reference"]] <- list(thresholds = refThresh, alarms = alarmRef)

# Extract thresholds and alarms for each method
thresholds <- sapply(results, "[[", "thresholds")
alarms <- sapply(results, "[[", "alarms")

# Analyse basic statistics of alarms
alarms.n <- sapply(alarms, length)
alarms.mean <- sapply(alarms, mean)

# Compare to a predefined cutpoints
cutPts <- seq(30, 60, 5)
sensitivity <- falseAlarms <- specificity <- youdenJ <- 
  precision <- F2 <- list()
for (i in 1:length(cutPts)){
  trueAlarms <- episodes(datatab$OM, cutPts[i])
  found <- lapply(alarms, function(x) trueAlarms[trueAlarms$t %in% names(x),])
  sensitivity[[i]] <- sapply(found, function(x) nrow(x) / nrow(trueAlarms))
  falseAlarms[[i]] <- sapply(alarms, 
    function(x) sum(!names(x) %in% trueAlarms$t))
  specificity[[i]] <- 1 - falseAlarms[[i]] / (n - nrow(trueAlarms))
  youdenJ[[i]] <- sensitivity[[i]] + specificity[[i]] - 1
  precision[[i]] <- sapply(alarms, function(x) 
    mean(as.integer(names(x)) %in% trueAlarms$t))
  F2[[i]] <- sapply(alarms, function(x) 
    Fscore(as.integer(names(x)), trueAlarms$t, 2))
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
bootRes[["CART"]] <- lapply(bsamples, function(dat) 
  CART.apply(dat$OM, dat[,indicator.names]))
bootRes[["MARS"]] <- lapply(bsamples, function(dat) 
  MARS.apply(dat$OM, dat[,indicator.names]))
bootRes[["PRIM"]] <- lapply(bsamples, function(dat) 
  PRIM.apply(dat$OM, dat[,indicator.names]))

# Obtain thresholds
bootThresh <- lapply(bootRes, sapply, "[[", "thresholds")

# Bootstrap estimates of bias and standard error
bootMean <- sapply(bootThresh, apply, 1, mean, na.rm = T)
bias <- bootMean - thresholds[, names(bootRes)]

bootSE <- sapply(bootThresh, apply, 1, sd, na.rm = T)

#---------------------------------------------------
#                 Plots
#---------------------------------------------------
respath <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/ArticleV5"

colPal <- brewer.pal(9, "Greys")
cols <- c("forestgreen", "cornflowerblue", "firebrick", "black")

# Visual representation of thresholds
x11() 
par(mar = c(5,4,4,7)+.1)
plot(datatab[,indicator.names], col = colPal[cut(datatab$OM, 9)], 
  xlab = which.X[1], ylab = which.X[2], cex.axis = 1.2, 
  cex.lab = 1.3, pch = 16, cex = 0.8, 
  xlim = c(10, max(datatab[,indicator.names[1]])),
  ylim = c(20, max(datatab[,indicator.names[2]])))
for(i in 1:length(results)){
  rect(max(thresholds[1,i], par("usr")[1] - 1, na.rm = T), 
    max(thresholds[2,i], par("usr")[3] - 1, na.rm = T), 
    par("usr")[2] + 1, par("usr")[4] + 1,
    border = cols[i], lwd = 2, lty = i)
}
image.scale(range(datatab$OM), colPal, zlab = "OM", 
  formatC.pars = list(digits = 0, format = "f"))
outerLegend("topcenter", names(results), lwd = 3, 
  lty = 1:(length(results) + 1), col = cols, ncol = 3, bty = "n")
  
dev.print(png, filename = sprintf("%s/Thresholds.png", respath), units = "in",
  res = 100)

# Sensitivity and specificity  
x11()
par(mfrow = c(3, 1), mar = c(5, 5, 4, 2) + .1)
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
matplot(cutPts, t(as.data.frame(F2)), type = "b", col = cols,
  pch = 14 + 1:length(results), xlab = "Cut point", ylab = expression(F[2]), 
  cex = 1.5, lwd = 3, cex.lab = 1.5, cex.axis = 1.2)
text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "c)", 
  cex = 3, xpd = T, adj = c(-0.5,1.2))
  
dev.print(png, filename = sprintf("%s/Performances.png", respath), units = "in",
  res = 100)

# Bootstrap results
x11()
par(mfrow = c(p, 1), cex.axis = 1.2, cex.lab = 1.3, mar = c(4, 5, 3, 2))
#layout(matrix(1:4, 4, 1, byrow = T), height = c(1, 4, 1, 4))
for (i in 1:p){    
  ylim <- range(c(lapply(bootThresh, "[", i, ), 
    bootMean[i,] + bootSE[i,], bootMean[i,] - bootSE[i,]),
    na.rm = T)
  vp <- vioplotMarc(lapply(bootThresh, "[", i, ), col = NA, border = cols, 
    pchMed = NA, pchMean = 16, rectBorder = NA, pchcex = 1,
    colMean = "black", pchMod = NA, drawMod = F, 
    lwd =2, lwdMed = NA, names = names(bootRes), ylim = ylim)
  title(ylab = "Threshold (°C)")
  abline(v = 2:length(bootRes) - 0.5, lty = 3, col = "darkgrey")
  abline(v = 1:length(bootRes))
  segments(1:length(bootRes) - 0.48, thresholds[i, names(bootRes)], 
    1:length(bootRes) + 0.48, thresholds[i, names(bootRes)],
    lwd = 3, lty = 2)
  arrows(1:length(bootRes), bootMean[i,], 1:length(bootRes),
    bootMean[i,] + bootSE[i,], angle = 90, length = .05,
    lwd = 2)
  arrows(1:length(bootRes), bootMean[i,], 1:length(bootRes),
    bootMean[i,] - bootSE[i,], angle = 90, length = .05,
    lwd = 2)
  title(main = sprintf("%s) %s", letters[i], which.X[i]))
}

dev.print(png, filename = sprintf("%s/Bootstrap_result.png", respath), units = "in",
  res = 100)