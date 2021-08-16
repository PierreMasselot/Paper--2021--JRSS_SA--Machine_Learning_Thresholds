################################################################################
#                              R code 
#                                                  
#     Machine learning approaches to identify thresholds in a heat-health 
#                      warning system context
#         Journal of the Royal Statistical Society - Series A
#                               2021
#
#                        Simulation functions
#
#                    Code Author: Pierre Masselot
#
################################################################################

# This script plots results from the simulation study. Must be executed after
#   11_Simulations.R

load("Results/11_Results_simulations.RData")

#----------------------------------------------------------------
#                         Useful objects
#----------------------------------------------------------------

# Number of models
nm <- ncol(sensitivity[[1]])

# Grid length for each parameter for plots
nr <- length(varParams$rho)
nt <- length(varParams$extType)

# Prepare X axis of plot
xvals <- 100 * unlist(varParams$extBetas)
nx <- length(xvals)
atx <- 1:(nm * nx) + rep(1:nx, each = nm)
atticks <- tapply(atx, rep(1:nx, each = nm), max)[-nx] + 1
atlabels <- tapply(atx, rep(1:nx, each = nm), mean)

# Model specific aesthetics (MOB, MARS, PRIM, AIM, GAM, SEG)
pal <- c("forestgreen", "cornflowerblue", "firebrick", "goldenrod", 
  "lightskyblue", "slategrey")
pchs <- 15:20

#----------------------------------------------------------------
#     Figure 1: F-score plot
#----------------------------------------------------------------

# Compute mean and 95 percent CI from simulated scores
medianf1 <- sapply(F1, apply, 2, mean, na.rm = T) 
quantlowf1 <- sapply(F1, apply, 2, quantile, .025, na.rm = T) 
quanthighf1 <- sapply(F1, apply, 2, quantile, .975, na.rm = T)

# Initialize plot
x11(width = 15)
par(oma = c(0, 0, 2, 0), mar = c(3, 4, 4, 2) + .1)
layout(rbind(matrix(1:(nr * nt), nrow = nr, ncol = nt), (nr*nt) + 1),
  height = c(.45, .45, .1))
for (i in seq_along(varParams$extType)){
  for (j in seq_along(varParams$rho)){
    # Select the experiments
    inds <- combParam[,"rho"] == varParams$rho[j] & 
      sapply(combParam[,"extType"], 
        function(x) identical(x, varParams$extType[[i]]))
    if1 <- c(medianf1[,inds])
    ilow <- c(quantlowf1[,inds])
    ihigh <- c(quanthighf1[,inds])
    
    # Plot results
    plot(atx, if1, pch = pchs, col = pal, cex = 1.5, 
      xlab = ifelse(j == nr, "Percent increase", ""), ylab = "F-score", 
      main = bquote(rho ~ "=" ~ .(varParams$rho[j])),
      ylim = c(0, 1), xaxt = "n")
    segments(x0 = atx, y0 = ilow, y1 = ihigh, lwd = 2, col = pal)
    abline(v = atticks, lty = 2)
    axis(1, at = atticks, labels = F)
    axis(1, at = atlabels, labels = xvals, lwd = 0)
    
    # Add labels
    if (j == nr) title(xlab = "Percent increase", xpd = NA) 
    if (j == 1) {
      mtext(sprintf("p = %i", i), side = 3, line = 4, 
        at = mean(par("usr")[1:2]), cex = 1.5, xpd = T)
    }
  }
}
# Add legend
par(mar = rep(0, 4))
plot.new()
legend("center", legend = colnames(sensitivity[[1]]), pch = pchs, lwd = 2,
  col = pal, ncol = nm, bty = "n", cex = 1.5)

# Save
dev.print(pdf, file = "Results/Figure1.pdf")

#----------------------------------------------------------------
#   Figure S1: Sensitivity
#----------------------------------------------------------------

# Compute mean and 95 percent CI from simulated scores
mediansens <- sapply(sensitivity, apply, 2, mean, na.rm = T) 
quantlowsens <- sapply(sensitivity, apply, 2, quantile, .025, na.rm = T) 
quanthighsens <- sapply(sensitivity, apply, 2, quantile, .975, na.rm = T)

# Initialize plot
x11(width = 15)
par(oma = c(0, 0, 2, 0), mar = c(3, 4, 4, 2) + .1)
layout(rbind(matrix(1:(nr * nt), nrow = nr, ncol = nt), (nr*nt) + 1),
  height = c(.45, .45, .1))
for (i in seq_along(varParams$extType)){
  for (j in seq_along(varParams$rho)){
    # Select the experiments
    inds <- combParam[,"rho"] == varParams$rho[j] & 
      sapply(combParam[,"extType"], 
        function(x) identical(x, varParams$extType[[i]]))
    isens <- c(mediansens[,inds])
    ilow <- c(quantlowsens[,inds])
    ihigh <- c(quanthighsens[,inds])
    
    # Plot results
    plot(atx, isens, pch = pchs, col = pal, cex = 1.5, 
      xlab = ifelse(j == nr, "Percent increase", ""), ylab = "Sensitivity", 
      main = bquote(rho ~ "=" ~ .(varParams$rho[j])),
      ylim = c(0, 1), xaxt = "n")
    segments(x0 = atx, y0 = ilow, y1 = ihigh, lwd = 2, col = pal)
    abline(v = atticks, lty = 2)
    axis(1, at = atticks, labels = F)
    axis(1, at = atlabels, labels = xvals, lwd = 0)
    
    # Add labels
    if (j == nr) title(xlab = "Percent increase", xpd = NA) 
    if (j == 1) {
      mtext(sprintf("p = %i", i), side = 3, line = 4, 
        at = mean(par("usr")[1:2]), cex = 1.5, xpd = T)
    }
  }
}
# Add legend
par(mar = rep(0, 4))
plot.new()
legend("center", legend = colnames(sensitivity[[1]]), pch = pchs, lwd = 2,
  col = pal, ncol = nm, bty = "n", cex = 1.5)

# Save
dev.print(pdf, file = "Results/FigSup1.pdf")

#----------------------------------------------------------------
#   Figure S2: Precision
#----------------------------------------------------------------

# Compute mean and 95 percent CI from simulated scores
medianprec <- sapply(precision, apply, 2, mean, na.rm = T) 
quantlowprec <- sapply(precision, apply, 2, quantile, .025, na.rm = T) 
quanthighprec <- sapply(precision, apply, 2, quantile, .975, na.rm = T)

# Initialize plot
x11(width = 15)
par(oma = c(0, 0, 2, 0), mar = c(3, 4, 4, 2) + .1)
layout(rbind(matrix(1:(nr * nt), nrow = nr, ncol = nt), (nr*nt) + 1),
  height = c(.45, .45, .1))
for (i in seq_along(varParams$extType)){
  for (j in seq_along(varParams$rho)){
    # Select the experiments
    inds <- combParam[,"rho"] == varParams$rho[j] & 
      sapply(combParam[,"extType"], 
        function(x) identical(x, varParams$extType[[i]]))
    iprec <- c(medianprec[,inds])
    ilow <- c(quantlowprec[,inds])
    ihigh <- c(quanthighprec[,inds])
    
    # Plot results
    plot(atx, iprec, pch = pchs, col = pal, cex = 1.5, 
      xlab = ifelse(j == nr, "Percent increase", ""), ylab = "Precision", 
      main = bquote(rho ~ "=" ~ .(varParams$rho[j])),
      ylim = c(0, 1), xaxt = "n")
    segments(x0 = atx, y0 = ilow, y1 = ihigh, lwd = 2, col = pal)
    abline(v = atticks, lty = 2)
    axis(1, at = atticks, labels = F)
    axis(1, at = atlabels, labels = xvals, lwd = 0)
    
    # Add labels
    if (j == nr) title(xlab = "Percent increase", xpd = NA) 
    if (j == 1) {
      mtext(sprintf("p = %i", i), side = 3, line = 4, 
        at = mean(par("usr")[1:2]), cex = 1.5, xpd = T)
    }
  }
}
# Add legend
par(mar = rep(0, 4))
plot.new()
legend("center", legend = colnames(sensitivity[[1]]), pch = pchs, lwd = 2,
  col = pal, ncol = nm, bty = "n", cex = 1.5)

# Save
dev.print(pdf, file = "Results/FigSup2.pdf")
