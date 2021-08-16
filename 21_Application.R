################################################################################
#                              R code 
#                                                  
#     Machine learning approaches to identify thresholds in a heat-health 
#                      warning system context
#         Journal of the Royal Statistical Society - Series A
#                               2021
#
#                        Main application
#
#                    Code Author: Pierre Masselot
#
################################################################################

# This script generates results shown in section 5.2. It details how to 
#   implement each considered method on real-world data, and computes basic
#   scores for each of them.

#----- Packages
library(earth) # MARS
library(mgcv) # GAM
library(partykit) # MOB
library(AIM) # AIM
library(primr) # PRIM (custom package, must be installed from github 
                # devtools::install_github("PierreMasselot/primr"))
library(splines) # Spline functions
library(colorspace) # Color palettes
library(hhws) # Functions for HHWS, must be installed from github 
              # devtools::install_github("PierreMasselot/hhws")
library(parallel) # Parallel computing
library(tools) # For compacting pdfs



library(forecast) # Time series functions

library(quantmod) 
library(RColorBrewer)
library(sm)
library(devtools)


      
#----- Load other functions
# Miscellaneous functions
source("00_Misc_functions.R")
# Wrappers to extract thresholds through the proposed methods
source("01_Threshold_functions.R")

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
cols <- c("forestgreen", "cornflowerblue", "firebrick", "goldenrod", "black")
pchs <- 15:18

# Cut points for OM
cutPts <- seq(30, 50, 5)

#---------------------------------------------------
#                 Data Loading
#---------------------------------------------------

# Load datafile
dataread <- read.table("Data.csv", header = T, sep=";")
dataread$Date <- as.POSIXlt(apply(dataread[,1:3], 1, paste, collapse = "-"))

# Extract summer months
datatab <- subset(dataread, Month %in% which.months)

# Create day-of-season variable
datatab$dos <- unlist(tapply(datatab[,1], datatab$Year, seq_along))

# Select variables
datatab <- datatab[,c(which.Y, which.X, "Year", "dos")]
datatab$Tmoy <- rowMeans(datatab[,which.X])

# Number of years, and days of data, and number of weather variables
nyear <- length(unique(datatab$Year))
n <- nrow(datatab)
p <- length(which.X)

#---------------------------------------------------
#               Prepare indicators
#---------------------------------------------------

L <- length(Lweights)
indicator.names <- which.X

# Compute the moving-average
indicators <- apply(datatab[,which.X], 2, stats::filter, 
  Lweights, sides = 1)
# To account for breaks in the series
indicators[which(diff(datatab$Date, L) > L) + L,] <- NA 

# Adding indicators to the data table
datatab[,indicator.names] <- indicators

# Remove NAs from data
datatab <- na.omit(datatab)
n <- nrow(datatab)

#---------------------------------------------------
#            Algorithms application
#---------------------------------------------------

# Initialize results
results <- list()

#----- Model-based regression trees (MOB) -----

# Grow the tree 
treefit <- glmtree(Death ~ Tmoy + ns(dos, 4) + ns(Year, round(nyear / 10)) | 
  Tmin + Tmax, data = datatab, family = "quasipoisson", minsize = 10)

# Extract thresholds from the tree object
results[["MOB"]] <- extractThresholds_party(treefit)$thresholds

# Infer best surrogate split for Tmin
surrfit <- glmtree(factor(Tmax > results[["MOB"]][2]) ~ Tmin, data = datatab, 
  family = "binomial", maxdepth = 2, minsize = 10)

# Plot the tree
x11(width = 17, height = 7)
plot(treefit, ip_args = list(id = FALSE, pval = FALSE), 
  terminal_panel = node_bivplotmod, ep_args = list(digits = 1), 
  tp_args = list(which = 1, lwd = 3, id = FALSE, linecol = cols[1], 
    pch = 16, pointcol = grey(.8)))

dev.print(pdf, file = "Results/Figure2.pdf")

#----- Multivariate adaptive regression splines (MARS) -----

# Fit MARS
marsfit <- earth(Death ~ Tmin + Tmax + dos + Year, 
  data = datatab, degree = p, endspan = 10, 
  glm = list(family = "quasipoisson"))

# Extract thresholds
results[["MARS"]] <- apply(marsfit$cuts[,indicator.names], 2, max)

## Plot the surface

# Coordinate grid for the surface
surfgrid <- do.call(expand.grid, 
  lapply(datatab[,indicator.names], 
    function(x) seq(min(x), max(x), length.out = 1000) )
)
datpred <- cbind(surfgrid, dos = mean(datatab$dos), Year = mean(datatab$Year))
surf <- matrix(predict(marsfit, datpred, type = "response"), nrow = 1000)

# Level plot
x11()
filled.contour(unique(surfgrid[,1]), unique(surfgrid[,2]), 
  surf, zlim = range(surf), 
  color.palette = function(n) sequential_hcl(n, rev = T), cex.axis = 1.2,
  key.title = title(main = "Death", cex = 1.3), 
  plot.title = {
    title(xlab = indicator.names[1], cex.lab = 1.3)
    title(ylab = indicator.names[2], cex.lab = 1.3)
    points(Tmax ~ Tmin, data = datatab, pch = 16, col = grey(.5))
    abline(v = unique(marsfit$cuts[marsfit$cuts[,1] != 0,1]), lty = 2)
    abline(v = results[["MARS"]][1], lwd = 2)
    mtext(formatC(results[["MARS"]][1], format = "f", digits = 1), 
      at = results[["MARS"]][1], line = 0.5, cex = 1.5)
    abline(h = unique(marsfit$cuts[marsfit$cuts[,2] != 0,2]), lty = 2)
    abline(h = results[["MARS"]][2], lwd = 2)
    mtext(formatC(results[["MARS"]][2], format = "f", digits = 1), 
      at = results[["MARS"]][2], line = 0.5, cex = 1.5, side = 2)
  },
)

# dev.print(pdf, file = "Results/Figure3.pdf")
dev.print(png, filename = "Results/Figure3.png", units = "in", res = 600)


#----- Patient rule-induction method (PRIM) -----

# Initial peeling
peelres <- peeling(datatab$Death, datatab[,indicator.names], 
  beta.stop = 10/n, peeling.side = -1,
  obj.fun = function(y, x, inbox){
    y <- y[inbox]
    x <- x[inbox,]
    dat <- data.frame(y, x, datatab[inbox, c("Tmoy", "dos", "Year")])
    fit <- glm(y ~ Tmoy + ns(dos, 4)+ ns(Year, round(nyear / 10)), 
      data = dat, family = "quasipoisson")
    pred <- coef(fit)[2]
    return(exp(pred))
})

# Choose number of observations above threshold (by analyzing peeling trajectory)
chosen <- 11

# Thresholds refinement by pasting 
primres <- pasting(peelres, support = chosen/n)

# Thresholds extraction
results[["PRIM"]] <- sapply(extract.box(primres)$limits[[1]], "[", 1) 

# Plotting the peeling trajectory
x11(width = 11)
plot_trajectory(peelres, pch = 16, col = cols[3], ylab = "RR increase rate", 
  xlim = c(0.002, 0.02), support = chosen/n, abline.pars = list(lwd = 2, lty = 2),
  ytype = "rel.diff", cex = 1.5)
mtext(sprintf("n = %i", chosen), at = chosen/n, cex = 1.5, line = 0.5)

dev.print(pdf, file = "Results/Figure4.pdf")

#----- Adaptive index models (AIM) -----

# Cross-validation to obtain the optimal number of steps
cvb <- cv.lm.main(datatab[,c("Tmin", "Tmax", "dos", "Year")], datatab$Death, 
  nsteps = 3 * 4, maxnumcut = 3, backfit = T, mincut = 10/n)

# Fit the AIM with the optimal number of steps
aimfit <- lm.main(datatab[,c("Tmin", "Tmax", "dos", "Year")], datatab$Death, 
  nsteps = cvb$kmax, maxnumcut = 3, backfit = T, mincut = 10/n)

# Extract thresholds
rules <- aimfit$res[[cvb$kmax]]
maxrules <- aggregate(rules[,2], by = list(var = rules[,1]), max)
results[["AIM"]] <- rep(NA_real_, p)
results[["AIM"]][maxrules$var[1:2]] <- maxrules$x[1:2]

# Predict index on the grid defined above for MARS surface
surfaim <- matrix(index.prediction(rules, datpred), nrow = 1000)
predscores <- sort(unique(c(surfaim)))

# Plot the surface
x11()
filled.contour(unique(surfgrid[,1]), unique(surfgrid[,2]), surfaim,
  color.palette = function(n) hcl.colors(n, "YlOrBr", rev = T), cex.axis = 1.2,
  key.title = title(main = "Score", cex = 1.3), levels = predscores,
  key.axes = axis(4, at = predscores + 0.5, labels = seq_along(predscores)),
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

# dev.print(pdf, file = "Results/Figure5.pdf")
dev.print(png, filename = "Results/Figure5.png", units = "in", res = 600)


#----- Saving estimated thresholds
thresholds <- Reduce(rbind, results)
rownames(thresholds) <- names(results)
write.table(thresholds, 
  file = "Results/Results_application.csv",
  quote = FALSE, sep = ";")


#---------------------------------------------------
#            Results analysis
#---------------------------------------------------

# Add reference (Chebana et al. 2013)
results[["Reference"]] <- c(20, 33)
  
# Compute over-mortality to define episodes
em.form <- sprintf("%s ~ ns(dos, 4) + ns(Year, round(nyear / 10))",
  which.Y)
EM <- lm(as.formula(em.form), datatab)$fitted.values
OM <- excess(datatab[,which.Y], EM) 

# Extract alarms for each method
alarms <- lapply(results, extract_alarms, x = datatab[,indicator.names], 
  y = OM) 

# Analyse basic statistics of alarms
alarms.n <- sapply(alarms, length)
episodes.n <- sapply(alarms, 
  function(x) sum(diff(as.integer(names(x))) > 3) + 1)
alarms.mean <- sapply(alarms, mean)
alarms.coverage <- alarms.mean * alarms.n / n

# Compare to a predefined cutpoint
sensitivity <- falseAlarms <- specificity <- precision <- F1 <- F2 <- list()
for (i in 1:length(cutPts)){
  trueAlarms <- episodes(OM, cutPts[i])
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
