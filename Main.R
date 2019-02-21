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

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/1 - Thresholds/Paper--Machine-Learning-Thresholds") # Data path
              
source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/PRIM/R/PRIM_functions.R")
source("Functions_rpartExtractResults.R")

#---------------------------------------------------
#                 Parameters
#---------------------------------------------------

# Parameters on data to consider
which.Y <- "Death"             # The response
which.X <- c("RH", "Tmax")   # The indicators
which.months <- 5:9            # the considered months

# Model parameters
L <- 2         # The lag for indicators

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

# Compute the moving-average
indicators <- apply(datatab[,which.X], 2, filter, 
  rep(1 / (L+1), L + 1), sides = 1)
# To account for breaks in the series
indicators[which(diff(datatab$Date, L) > L) + L,] <- NA 

# Adding indicators to the data table
datatab[,sprintf("Z%s", which.X)] <- indicators

#---------------------------------------------------
#            Remove NAs from data table
#---------------------------------------------------

datatab <- na.omit(datatab)

####################################################
#                   CART
####################################################

library(rpart) # For CART
library(rpart.plot)  # For plotting the result

source("Functions_rpartExtractResults.R")

#---------------------------------------------------
#                CART application
#---------------------------------------------------

# Tree growing
# The argument cp = 0.0001 ensures that a large tree is grown  
# (Note: rpart randomly select a subset of observations to grow the tree for performances purposes. The result may thus be slightly different than shown in the manuscript).
tree.growing <- rpart(OM ~ ., data = data.frame(OM = OM, indicators), method = "anova", control = rpart.control(minsplit = 10, cp = 0.0001))

# Determine the best subtree (through its complexity parameter 'cp')
minind <- which.min(tree.growing$cptable[,"xerror"])
best_cp <- tree.growing$cptable[minind, "CP"]

# Prune the tree
final.tree <- prune(tree.growing, best_cp)

#--------------------------------------------------
#                   Results
#--------------------------------------------------

# Obtain the thresholds
CART.result <- get.box(final.tree, varnames = c("Tmin", "Tmax"))
CART.thresholds <- CART.result$box[1,]
# Threshold for Tmin is equal in this case to its minimum value, which is equivalent to no threshold

# Alerts
CART.alerts <- indicators[,1] >= CART.result$box[1,1] & indicators[,2] >= CART.result$box[1,2]

#--------------------------------------------------
#                   Plots
#--------------------------------------------------

## Tree

cols <- rep("black", nrow(final.tree$frame))
cexs <- rep(0.8, nrow(final.tree$frame))
lwds <- rep(1, nrow(final.tree$frame))

path <- colpath_to_node(CART.result$node, final.tree)
cols[path] <- "forestgreen"
cexs[path] <- 1.5
lwds[path] <- 4 

x11() # Figure 2a)
prp(final.tree, extra = 1, left = FALSE, fallen.leaves = F, col=cols, branch.col=cols, split.col=cols, cex = cexs, branch.lwd = lwds, digits = 3)

## Partition

x11() # Figure 2b)
par(mar = c(5,4,4,4)+.1)
plot(indicators, col = "white", xlab = "Tmin", ylab = "Tmax", cex.axis = 1.2, cex.lab = 1.3, pch = 16, cex = 0.8) # Initialize an empty plot
rect(CART.result$box[1,1], CART.result$box[1,2], par("usr")[2]+1, par("usr")[4]+1, col = "forestgreen", border = NA) # Add background for the alert subset
points(indicators, pch = 16, cex = 0.8, col = ifelse(CART.alerts, "black", "darkgrey")) # Add datapoints
add.partition(final.tree, segment.pars = list(lty = 2, col = "black", lwd = 2), leaf.text = "none")
mtext(side = 4, at = CART.result$box[1,2], text = round(CART.result$box[1,2], digits = 1), col = "forestgreen", cex = 1.5, line = .5, las = 1) # Add text for threshold on Tmax
box()


####################################################
#                   MARS
#################################################### 

library(earth) # contains MARS

#---------------------------------------------------
#                MARS application
#---------------------------------------------------

MARS.result <- earth(x = indicators, y = OM)

#---------------------------------------------------
#                Extracting the result
#---------------------------------------------------

#Threshold
MARS.thresholds <- apply(MARS.result$cuts, 2, max) # Taking the two extremal cuts

# Alerts
MARS.alerts <- indicators[,1] >= MARS.thresholds[1] & indicators[,2] >= MARS.thresholds[2]

#--------------------------------------------------
#                     Plot
#--------------------------------------------------

x11() # Figure 3
par(mar = c(5,4,4,4)+.1)
plot(indicators, col = "white", xlab = "Tmin", ylab = "Tmax", cex.axis = 1.2, cex.lab = 1.3, pch = 16, cex = 0.8) # Initialize an empty plot
rect(MARS.thresholds[1], MARS.thresholds[2], par("usr")[2]+1, par("usr")[4]+1, col = "cornflowerblue", border = NA) # Add background for the alert subset
points(indicators, pch = 16, cex = 0.8, col = ifelse(MARS.alerts, "black", "darkgrey"))
abline(v = unique(MARS.result$cuts[MARS.result$cuts[,1]!=0,1]), lty = 2, lwd = 2) # Add cuts on Tmin
abline(h = unique(MARS.result$cuts[MARS.result$cuts[,2]!=0,2]), lty = 2, lwd = 2) # Add cuts on Tmax
mtext(side = 3, at = MARS.thresholds[1], text = round(MARS.thresholds[1], digits = 1), col = "cornflowerblue", cex = 1.5, line = .5) # Add text for threshold on Tmin
mtext(side = 4, at = MARS.thresholds[2], text = round(MARS.thresholds[2], digits = 1), col = "cornflowerblue", cex = 1.5, line = .5, las = 1) # Add text for threshold on Tmax
box()


####################################################
#                   PRIM
#################################################### 

source("Functions_PRIM.R")

#---------------------------------------------------
#                PRIM application
#---------------------------------------------------

# Apply the peeling algorithm with a peeling fraction of 5%, a stopping criterion of 6 observations and the constraint that only the box can shrink only through its lower sides.
peeling.result <- peeling.sequence(OM, indicators, alpha = .05, beta.stop = 6/n, peeling.side = "left")

# Extract the final box
chosen.support <- 10 # In this case, the box containing 10 observations is chosen since the average OM only slightly increases in further steps
peeled.box <- extract.box(peeling.result, chosen.support/n)

# Optional: refine the boundaries of the box by pasting
PRIM.result <- pasting.sequence(OM, indicators, small.box = peeled.box$limits, peeling.side = "left", alpha = 0.01)

#---------------------------------------------------
#                Extracting the results
#---------------------------------------------------

#Thresholds
PRIM.thresholds <- PRIM.result$limits[1,]

# Alerts
PRIM.alerts <- indicators[,1] >= PRIM.thresholds[1] & indicators[,2] >= PRIM.thresholds[2]

#--------------------------------------------------
#                     Plots
#--------------------------------------------------

# Plot of the peeling sequence to choose the best box
x11() # Figure 4a)
plot(round(peeling.result$support*n), peeling.result$yfun, type = "b", xlab = "Observation number", ylab = "Average Over-mortality", pch = 21, bg = c("white","firebrick")[1+(peeling.result$support*n <= chosen.support)], col = "black", log = "x", cex.axis = 1.2, cex.lab = 1.3, xlim = rev(range(peeling.result$support*n)))
abline(v = chosen.support, lty = 2, lwd = 3, col = "firebrick")
mtext(side = 3, at = chosen.support, text = chosen.support, cex = 1.5, col = "firebrick", line = .5)

# Plot of the box
x11() # Figure 4b)
par(mar = c(5,4,4,4)+.1)
plot(indicators, col = "white", xlab = "Tmin", ylab = "Tmax", cex.axis = 1.2, cex.lab = 1.3, pch = 16, cex = 0.8) # Initialize an empty plot
rect(PRIM.thresholds[1], PRIM.thresholds[2], par("usr")[2]+1, par("usr")[4]+1, col = "firebrick", lty = 2, lwd = 2) # Add background for the alert subset
points(indicators, pch = 16, cex = 0.8, col = ifelse(PRIM.alerts, "black", "darkgrey"))
mtext(side = 3, at = PRIM.thresholds[1], text = round(PRIM.thresholds[1], digits = 1), col = "firebrick", cex = 1.5, line = .5) # Add text for threshold on Tmin
mtext(side = 4, at = PRIM.thresholds[2], text = round(PRIM.thresholds[2], digits = 1), col = "firebrick", cex = 1.5, line = .5, las = 1) # Add text for threshold on Tmax
box()

####################################################
#       Comparison with classical method
#################################################### 

library(Kendall)

source("Functions_classicalMethod.R")

#--------------------------------------------------
#        Determine over-mortality episodes
#--------------------------------------------------

# Estimation of OM trend
MK.OM <- MannKendall(OM) # The trend is not significant
trend.OM <- coef(lm(OM ~ t, data = data.frame(OM = OM, t = 1:n)))[2] # Estimation of the trend. Not used in the following

# extract excesses
OMT <- episodes(OM, u = 40, covariates = X[tt,], uc = c(0, 28), l = 3) # Over-mortality threshold fixed at 40%


#--------------------------------------------------
#             Find thresholds
#--------------------------------------------------

# Try different thresholds combinations. We force the weights to be 1/3 here for consistency with other methods
classical.result <- find.threshold(aperm(tmp[tt,,], c(1,3,2)), episodes=OMT$excesses, alphas = 1/3, u.grid=list(10:20,30:40), trim = 30)

# Choose best thresholds according to sensitivity and specificity
chosen.index <- 9
classical.thresholds <- as.numeric(classical.result[9, colnames(X)])

# Alerts
classical.alerts <- indicators[,1] >= classical.thresholds[1] & indicators[,2] >= classical.thresholds[2]


#--------------------------------------------------
#        Comparison with proposed methods
#--------------------------------------------------

# gather all alerts into one list
alerts.list <- list(CART = OM[CART.alerts], MARS = OM[MARS.alerts], PRIM = OM[PRIM.alerts], Classical = OM[classical.alerts])

# Boxplot of the alerts OM
x11() #Figure5a)
boxplot(alerts.list, border = c("forestgreen", "cornflowerblue", "firebrick", grey(.3)), lwd = 2, ylim = range(OM), ylab = "Over-mortality", cex.lab = 1.3, cex.axis = 1.2, xlab = "", varwidth = T)


####################################################
#       Calibration/validation design
#################################################### 

# Inspired by Hajat, S. et al. Heat–Health Warning Systems: A Comparison of the Predictive Capacity of Different Approaches to Identifying Dangerously Hot Days. American journal of public health 100, 1137–1144 (2010).

#--------------------------------------------------
#      Separate calibration/validation samples
#--------------------------------------------------

# Order summers' temperatures
Tmean <- apply(indicators, 1, mean) # We use mean temperature to sort summers' heat
summer_heat <- aggregate(Tmean, by = list(year = dates$year), mean) # mean temperature of each year. 
summer_ord <- order(summer_heat[,2], decreasing = T) # decreasing order

# Attribute odd-ordered years to the calibration sample and even-ordered years to the validation sample
sample.ind <- rep_len(1:2, nyear)
names(sample.ind) <- summer_heat[summer_ord,1] # years with a 1 are attributed to the calibration sample

# Attribute each day to a sample
sample.day <- sample.ind[as.character(dates$year)] # attribute 1 or 2 to each day
# calibration
indicators.calib <- indicators[sample.day == 1,]
OM.calib <- OM[sample.day == 1]
# validation
indicators.valid <- indicators[sample.day == 2,]
OM.valid <- OM[sample.day == 2]

#--------------------------------------------------
#                      CART
#--------------------------------------------------

# calibrate CART
tree.calib <- rpart(OM ~ ., data = data.frame(OM = OM.calib, indicators.calib), method = "anova", control = rpart.control(minsplit = 10, cp = 0.0001))
minind <- which.min(tree.calib$cptable[,"xerror"])
best_cp <- tree.calib$cptable[minind, "CP"]
final.tree.calib <- prune(tree.calib, best_cp)
CART.calib <- get.box(final.tree.calib, varnames = c("Tmin", "Tmax"))
CART.thresholds.calib <- CART.calib$box[1,]

# Alerts on validation sample
CART.alerts.valid <- indicators.valid[,1] >= CART.thresholds.calib[1] & indicators.valid[,2] >= CART.thresholds.calib[2]

#--------------------------------------------------
#                      MARS
#--------------------------------------------------

# calibration
MARS.calib <- earth(x = indicators.calib, y = OM.calib)
MARS.thresholds.calib <- apply(MARS.calib$cuts, 2, max) # Taking the two extremal cuts

# Alerts
MARS.alerts.valid <- indicators.valid[,1] >= MARS.thresholds.calib[1] & indicators.valid[,2] >= MARS.thresholds.calib[2]

#--------------------------------------------------
#                      PRIM
#--------------------------------------------------

# Calibration
peeling.calib <- peeling.sequence(OM.calib, indicators.calib, alpha = .05, beta.stop = 6/n, peeling.side = "left")
peeled.calib <- extract.box(peeling.calib, 8/length(OM.calib)) # chosen after inspecting the peeling trajectory
PRIM.calib <- pasting.sequence(OM.calib, indicators.calib, small.box = peeled.calib$limits, peeling.side = "left", alpha = 0.01)
PRIM.thresholds.calib <- PRIM.calib$limits[1,]

# Validation
PRIM.alerts.valid <- indicators.valid[,1] >= PRIM.thresholds.calib[1] & indicators.valid[,2] >= PRIM.thresholds.calib[2]

#--------------------------------------------------
#                Classical method
#--------------------------------------------------

# Calibration
# We do not look for a trend since none was found in the whole sample
OMT.calib <- episodes(OM.calib, u = 30, covariates = X[tt,][sample.day == 1,], uc = c(0, 28), l = 3) # u chosen at 30 here to keep an adequate number of OM episodes
classical.calib <- find.threshold(aperm(tmp[tt,,][sample.day == 1,,], c(1,3,2)), episodes=OMT.calib$excesses, alphas = 1/3, u.grid=list(10:20,30:40), trim = 30)
classical.thresholds.calib <- as.numeric(classical.calib["12", colnames(X)])

# Validation
classical.alerts.valid <- indicators.valid[,1] >= classical.thresholds.calib[1] & indicators.valid[,2] >= classical.thresholds.calib[2]

#--------------------------------------------------
#                Boxplots
#--------------------------------------------------

valid.list <- list(CART = OM.valid[CART.alerts.valid], MARS = OM.valid[MARS.alerts.valid], PRIM = OM.valid[PRIM.alerts.valid], Classical = OM.valid[classical.alerts.valid])

x11() #Figure5b)
boxplot(valid.list, border = c("forestgreen", "cornflowerblue", "firebrick", grey(.3)), lwd = 2, ylim = range(OM), ylab = "Over-mortality", cex.lab = 1.3, cex.axis = 1.2, xlab = "", varwidth = T)