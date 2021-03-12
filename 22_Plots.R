#########################################################
#            Main code to produce results in
#                                                  
#     "Objective methods for threshold assessment in 
#          weather-health watch warning systems"
#
#             Author: Pierre Masselot
#
#                    Journal
#
#              Plots
#
#########################################################

load("Results/21_BootRes.RData")

#-----------------------------
# Thresholds
#-----------------------------

# Extract all thresholds
boot_thresh <- sapply(bootRes, sapply, "[[", "thresholds", simplify = "array")

# Object for plot
nm <- length(bootRes)
ylim <- range(boot_thresh, na.rm = TRUE)
xlim <- c(.5, p * nm + 0.5)

# Create device window
x11(width = 10, height = 6)
par(yaxs = "i", cex.axis = 1.2, cex.lab = 1.3, mar = c(5, 2, 4, 1) + .1)
layout(matrix(1:2, nrow = 1), width = c(0.7, 0.3))

# Threhsolds
for (j in 1:p){    
  resj <- boot_thresh[j,,]
  boxplot(resj, horizontal = T, at = (1:nm) * p - j + 1, yaxt = "n",
    border = cols, varwidth = T, col = NA,
    lwd = 2, ylim = ylim, xlim = xlim, add = j > 1, pch = j, 
    xlab = "Temperature", main = "a) Bootstrap thresholds")
  obsThresh <- sapply(results[1:nm], "[", j)
  segments(obsThresh, (1:nm - 1) * p + 0.6, 
    obsThresh, (1:nm) * p + 0.4, 
    lwd = 3, lty = 1, lend = 2)
}
abline(h = 0.5 + p * 1:(nm - 1), lty = 3)
axis(side = 2, at = (1:nm) * p - (p - 1)/2, labels = names(bootRes),
  tick = F)

# Number of unfound thresholds
percentNA <- apply(boot_thresh, c(1,3), function(x) mean(is.na(x)))
bp <- barplot(percentNA[p:1,], beside = TRUE, horiz = TRUE, xlim = c(0,1),
  col = rep(cols[1:nm], each = p), space = c(0, 0.1), xaxt = "n",
  xlab = "Proportion (%)", main = "b) Discarded variables")
axis(1, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2) * 100)
abline(h = bp[1, -1] - 0.55, lty = 3)
box()

# Save
dev.print(png, filename = "Results/Fig6_Boot_thresh.png", 
  units = "in", res = 100)
dev.copy2eps(file = "Results/Fig6_Boot_thresh.eps")

#-----------------------------
# Detection of over-mortality events
#-----------------------------

# Computation of scores for each method
boot_sensitivity <- boot_precision <- boot_F1 <- boot_F2 <- 
  vector("list", length(cutPts))
for (i in 1:length(cutPts)){
  # True extreme days
  trueX <- lapply(bsamples, function(b) which(OM[b] > cutPts[i]))
  ntrue <- sapply(trueX, length)
  
  # Compute scores for each simulation of each model
  scores <- lapply(bootRes, function(x){
    preddays <- lapply(x, "[[", "alerts")
    found <- mapply("%in%", preddays, trueX, SIMPLIFY = F)
    sens <- sapply(found, sum) / ntrue
    prec <- sapply(found, sum) / sapply(preddays, length)
    prec[prec == Inf] <- 0
    F1 <- mapply(Fscore, preddays, trueX, 1)
    F2 <- mapply(Fscore, preddays, trueX, 2)
    list(sensitivity = sens, precision = prec, F1 = F1, F2 = F2)
  })
  
  # Store results
  boot_sensitivity[[i]] <- sapply(scores, "[[", "sensitivity")
  boot_precision[[i]] <- sapply(scores, "[[", "precision")
  boot_F1[[i]] <- sapply(scores, "[[", "F1")
  boot_F2[[i]] <- sapply(scores, "[[", "F2")
}



#---- Plot results ----

# Params
quant <- c(0.025, 0.975)

# useful objects
nm <- length(bootRes)
nx <- length(cutPts)
atx <- 1:(nm * nx) + rep(1:nx, each = nm)
atticks <- tapply(atx, rep(1:nx, each = nm), max)[-nx] + 1
atlabels <- tapply(atx, rep(1:nx, each = nm), mean)

x11(height = 10)
layout(matrix(1:4, ncol = 1), height = c(rep(4, 3), 1))
par(mar = c(5, 5, 4, 2) + .1, cex.main = 2)

#--- Sensitivity
# Compute confidence intervals
sens_low <- sapply(boot_sensitivity, apply, 2, quantile, quant[1], na.rm = T)
sens_high <- sapply(boot_sensitivity, apply, 2, quantile, quant[2], na.rm = T) 
# Plot
plot(atx, sapply(sensitivity, "[", 1:4), pch = pchs, col = cols[1:4], cex = 1.5, 
  xlab = "", ylab = "Sensitivity", ylim = range(c(sens_low, sens_high)), 
  xaxt = "n", main = "a) Sensitivity")
segments(x0 = atx, y0 = sens_low, y1 = sens_high, lwd = 2, col = cols[1:4])
abline(v = atticks, lty = 2)
axis(1, at = atticks, labels = F)
axis(1, at = atlabels, labels = cutPts, lwd = 0)
segments(x0 = c(0, atticks), x1 = c(atticks, max(atx) + 1), 
  y0 = sapply(sensitivity, "[", 5), lend = 3, lwd = 2)
# text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "a) Sensitivity", 
#   cex = 3, xpd = T, adj = c(-0.5,1.2))

#--- Precision
# Compute confidence intervals
prec_low <- sapply(boot_precision, apply, 2, quantile, quant[1], na.rm = T)
prec_high <- sapply(boot_precision, apply, 2, quantile, quant[2], na.rm = T) 
# Plot
plot(atx, sapply(precision, "[", 1:4), pch = pchs, col = cols[1:4], cex = 1.5, 
  xlab = "", ylab = "Precision", ylim = range(c(prec_low, prec_high)), 
  xaxt = "n", main = "b) Precision")
segments(x0 = atx, y0 = prec_low, y1 = prec_high, lwd = 2, col = cols[1:4])
abline(v = atticks, lty = 2)
axis(1, at = atticks, labels = F)
axis(1, at = atlabels, labels = cutPts, lwd = 0)
segments(x0 = c(0, atticks), x1 = c(atticks, max(atx) + 1), 
  y0 = sapply(precision, "[", 5), lend = 3, lwd = 2)
# text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "b) Precision", 
#   cex = 3, xpd = T, adj = c(-0.5,1.2))

#--- F1
# Compute confidence intervals
F1_low <- sapply(boot_F1, apply, 2, quantile, quant[1], na.rm = T)
F1_high <- sapply(boot_F1, apply, 2, quantile, quant[2], na.rm = T) 
# Plot
plot(atx, sapply(F1, "[", 1:4), pch = pchs, col = cols[1:4], cex = 1.5, 
  xlab = "", ylab = "F1 score", ylim = range(c(F1_low, F1_high)), 
  xaxt = "n", main = "c) F1")
segments(x0 = atx, y0 = F1_low, y1 = F1_high, lwd = 2, col = cols[1:4])
abline(v = atticks, lty = 2)
axis(1, at = atticks, labels = F)
axis(1, at = atlabels, labels = cutPts, lwd = 0)
segments(x0 = c(0, atticks), x1 = c(atticks, max(atx) + 1), 
  y0 = sapply(F1, "[", 5), lend = 3, lwd = 2)
# text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "c) F1", 
#   cex = 3, xpd = T, adj = c(-0.5,1.2))

#--- F2
# Compute confidence intervals
# F2_low <- sapply(boot_F2, apply, 2, quantile, quant[1], na.rm = T)
# F2_high <- sapply(boot_F2, apply, 2, quantile, quant[2], na.rm = T) 
# # Plot
# plot(atx, sapply(F2, "[", 1:4), pch = pchs, col = cols[1:4], cex = 1.5, 
#   xlab = "", ylab = "Sensitivity", ylim = range(c(F2_low, F2_high)), 
#   xaxt = "n", main = "d) F2")
# segments(x0 = atx, y0 = F2_low, y1 = F2_high, lwd = 2, col = cols[1:4])
# abline(v = atticks, lty = 2)
# axis(1, at = atticks, labels = F)
# axis(1, at = atlabels, labels = cutPts, lwd = 0)
# title(xlab = "Cut point")
# segments(x0 = c(0, atticks), x1 = c(atticks, max(atx) + 1), 
#   y0 = sapply(F2, "[", 5), lend = 3, lwd = 2)
# text(line2user(par("mar")[2], 2), line2user(par("mar")[3], 3), "d) F2", 
#   cex = 3, xpd = T, adj = c(-0.5,1.2))

#--- Add legend
par(mar = rep(0, 4))
plot.new()
legend("center", legend = names(results), col = cols, pch = c(pchs, NA),
  lty = c(rep(NA, nm), 1), bty = "n", ncol = nm + 1, lwd = 2, cex = 1.5)

dev.print(png, filename = "Results/Fig7_Boot_scores.png", 
  units = "in", res = 100)
dev.copy2eps(file = "Results/Fig7_Boot_scores.eps")


#-----------------------------
# Table 2
#-----------------------------

# Thresholds
apply(boot_thresh, c(1,3), quantile, c(.025, .975), na.rm = T)

preddays <- lapply(bootRes, lapply, "[[", "alerts")
# Number of days
ndays <- sapply(preddays, sapply, length)
apply(ndays, 2, quantile, c(.025, .975))

# Number of episodes
episodes.n <- sapply(preddays, sapply, 
  function(x) sum(diff(x) > 3) + 1)
apply(episodes.n, 2, quantile, c(.025, .975))

# OM mean
alarms.mean <- matrix(NA, B, nm)
for (i in 1:nm){
  alarms.mean[,i] <- mapply(function(x, b) mean(OM[b][x]), 
    preddays[[i]], bsamples)
}
apply(alarms.mean, 2, quantile, c(.025, .975))
  
# Coverage
alarms.coverage <- 100 * mapply("*", as.data.frame(alarms.mean), as.data.frame(ndays), 
  SIMPLIFY = "matrix") / n
apply(alarms.coverage, 2, quantile, c(.025, .975))