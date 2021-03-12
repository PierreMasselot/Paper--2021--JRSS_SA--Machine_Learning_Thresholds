#########################################################
#
#            Parametric simulations in paper
#                 Plots of the results
#                                                  
#########################################################

load("Results/11_Results_simulations.RData")

#----------------------------------------------------------------
#                         Useful objects
#----------------------------------------------------------------

# Number of models
nm <- ncol(sensitivity[[1]])

# Number of each paremeters for plots
nr <- length(varParams$rho)
nt <- length(varParams$extType)

# X axis of plot
xvals <- 100 * unlist(varParams$extBetas)
nx <- length(xvals)
atx <- 1:(nm * nx) + rep(1:nx, each = nm)
atticks <- tapply(atx, rep(1:nx, each = nm), max)[-nx] + 1
atlabels <- tapply(atx, rep(1:nx, each = nm), mean)

# Model specific aesthetics
pal <- c("forestgreen", "cornflowerblue", "firebrick", "goldenrod", 
  "lightskyblue", "slategrey")
pchs <- 15:20

#----------------------------------------------------------------
#                         F1-score plot
#----------------------------------------------------------------

# Compute summary statistics of the criterion
medianf1 <- sapply(F1, apply, 2, mean, na.rm = T) 
quantlowf1 <- sapply(F1, apply, 2, quantile, .025, na.rm = T) 
quanthighf1 <- sapply(F1, apply, 2, quantile, .975, na.rm = T)

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
      xlab = ifelse(j == nr, "Percent increase", ""), ylab = expression(F[1]), 
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

dev.print(png, filename = "Results/Fig1_Simu_F1.png", 
  units = "in", res = 100)
dev.copy2eps(file = "Results/Fig1_Simu_F1.eps")

#----------------------------------------------------------------
#                         Sensitivity
#----------------------------------------------------------------

# Compute summary statistics of the criterion
mediansens <- sapply(sensitivity, apply, 2, mean, na.rm = T) 
quantlowsens <- sapply(sensitivity, apply, 2, quantile, .025, na.rm = T) 
quanthighsens <- sapply(sensitivity, apply, 2, quantile, .975, na.rm = T)

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
      xlab = ifelse(j == nr, "Percent increase", ""), ylab = expression(F[1]), 
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

dev.print(png, filename = "Results/FigSup1_Simu_sens.png", 
  units = "in", res = 100)

#----------------------------------------------------------------
#                         Sensitivity
#----------------------------------------------------------------

# Compute summary statistics of the criterion
medianprec <- sapply(precision, apply, 2, mean, na.rm = T) 
quantlowprec <- sapply(precision, apply, 2, quantile, .025, na.rm = T) 
quanthighprec <- sapply(precision, apply, 2, quantile, .975, na.rm = T)

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
      xlab = ifelse(j == nr, "Percent increase", ""), ylab = expression(F[1]), 
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

dev.print(png, filename = "Results/FigSup2_Simu_prec.png", 
  units = "in", res = 100)