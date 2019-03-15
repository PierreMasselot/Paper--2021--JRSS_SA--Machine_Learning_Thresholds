#########################################################
#
#            Parametric simulations in paper
#                                                  
#########################################################

library(quantmod)
library(rpart)
library(earth)
library(hhws)
library(lattice)

source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/PRIM/R/PRIM_functions.R")
source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/1 - Thresholds/Paper--Machine-Learning-Thresholds/Threshold_functions.R")
source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/1 - Thresholds/Paper--Machine-Learning-Thresholds/Functions_rpartExtractResults.R")
source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Export_functions.R")
source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Util_functions.R")
source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Plot_functions.R")

#---------------------------------------------------
#                 Relationship functions
#---------------------------------------------------

# Dose-response functions
fconstant <- function(x, Betas) Betas
flinear <- function(x, Betas) Betas[1] + Betas[2] * x
fJshape <- function(x, Betas){
  pow <- outer(x, 0:4, "^")
  pow %*% Betas
}

# Lag-response functions
wone <- function(l, Betas) ifelse(l == 0, 1, 0)
wconstant <- function(l, Betas) Betas
wdecay <- function(l, Betas) Betas[1] * exp(Betas[2] * l) 

#---------------------------------------------------
#              Data generating function
#---------------------------------------------------
#' Generate data for testing threshold finding methods
#'
#' @param B Number of simulated datasets.
#' @param n Record length of simulated datasets.
#' @param p Number of indicator variables.
#' @param Lsim Maximum lag for the relationship between response and indicators.
#' @param obetas Weights of each indicator's effect on the response.
#' @param ffuns Dose-response functions. Recycled if necessary.
#' @param fbetas Parameters of the dose-response functions.
#' @param wfuns Lag-response functions. Recycled if necessary.
#' @param wbetas Parameters or lag-response functions.
#' @param s Thresholds to be found. Must be between 0 and 1.
#' @param stype If "quantile" s corresponds to the quanile of X. If not, directly
#'   directly corresponds to the threshold.
#' @param extBetas Parameters for the linear function of extremes.
#' @param XdepOrder AR order for temporal dependence of indicators.
#' @param YdepOrder AR ordre for temporal dependence of response.
#' @param noise.sd Standard deviation of the random part of response, relative
#'    to the standard deviation of the deterministic part.
generate.data <- function(B = 1000, n = 5000, p = 2, Lsim = 1, 
  obetas = 1, ffuns = "fconstant", fbetas = list(), wfuns = "wone", 
  wbetas = list(), s = .8, stype = c("quantile", "absolute"), extBetas = 1,
  XdepOrder = 0, YdepOrder = 0, noise.sd = .2)
{
  params <- as.list(environment())
  
  # Recycling necessary parameters
  ffuns <- rep_len(ffuns, p)
  wfuns <- rep_len(wfuns, p)
  s <- rep_len(s, p)
  obetas <- rep_len(obetas, p + 1)
  wbetas <- rep_len(wbetas, p + 1)
  extBetas <- rep_len(extBetas, p + 1)
  
  # Predictors X
  suppressWarnings(Xsim <- replicate(p, 
    arima.sim(list(ar = XdepOrder), n, rnorm)))
  stype <- match.arg(stype)
  if (stype == "quantile"){
    Xsim <- apply(Xsim, 2, function(x) rank(x) / n)
  } else {
    Xsim <- apply(Xsim, 2, scale)
  }
  
  # Mean part
  Ypart <- Xsim
  fTrue <- overall <- list()
  for (j in 1:p){
    Xlags <- sapply(0:Lsim, Lag, x = Xsim[,j])
    Xfun <- function(x, l){do.call(ffuns[j], list(x, fbetas[[j]])) * 
      do.call(wfuns[j], list(l, wbetas[[j]]))
    }
    fTrue[[j]] <- outer(seq(0, 1, length.out = 20), 0:Lsim, Xfun)
    overall[[j]] <- apply(fTrue[[j]], 1, sum)
    Ypart[,j] <- apply(Xlags, 1, function(x) sum(Xfun(x, 0:Lsim)))
  } 
  Ypart <- apply(Ypart, 2, scale)
  
  linPred <- cbind(1, Ypart) %*% obetas
  linPred <- scale(linPred)
  
  # Extreme part
  uni.extremes <- mapply(">=", as.data.frame(Xsim), s)
  extremes <- apply(uni.extremes, 1, all)
  
  extBetas <- extBetas
  Xext <- mapply("-", as.data.frame(Xsim), s)
  extPred <- ifelse(extremes, cbind(1, Xext) %*% extBetas, 0)

  # Response simulation
  suppressWarnings(Ysim <- replicate(B, scale(linPred + extPred) + 
    arima.sim(list(ar = YdepOrder), n, rnorm, sd = noise.sd), simplify = TRUE))
  
  output <- list(Ysim = Ysim, Xsim = Xsim, linPred = linPred, extPred = extPred,
    fTrue = fTrue, overall = overall, parameters = params, s = s)
  return(output)
}


##############################################################################
#
#                            Simulations 

sim.name <- "BiLinear_XLinear_noDep"
SAVE <- FALSE
savepath <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/Simulations"

# Important constant parameters
p <- 2
s <- .95  # Quantile
B <- 1000

# Varying parameters between simulations
varParams <- within(list(), {
  Lsim <- 0:5; names(Lsim) <- sprintf("L%s", Lsim)
  extType <- list(jump = c(1, 0, 0), linear = c(0, 1, 1))
  extBetas <- list(X10 = .1, X30 = .3, X50 = .5, X100 = 1, X300 = 3, X500 = 5, 
                   X1000 = 10)
})

combParam <- do.call(expand.grid, varParams)
nc <- nrow(combParam)

# Other useful objects
cols <- c("forestgreen", "cornflowerblue", "firebrick")

# Scores of interest 
biasRes <- vector("list", nc)
SEres <- vector("list", nc)
RMSEres <- vector("list", nc)

#---- Loop for simulations ---
for(i in 1:nc){  
  print(sprintf("%i / %i", i, nc)); flush.console()
  
  #---- Generate data
  L <- combParam[i, "Lsim"]
  XBetas <- combParam[[i, "extType"]] * combParam[[i, "extBetas"]]
  Xdep <- combParam[i, "Xdep"]
  
  dataSim <- generate.data(B = B, Lsim = L, p = p, s = s,
    obetas = c(0, 1, 1), ffuns = "flinear", fbetas = list(c(0,1), c(0,1)), 
    wfuns = c("wconstant", "wdecay"), 
    wbetas = list(1 / (L + 1), c(1, -5 / (L + 1))),
    extBetas = XBetas, XdepOrder = .6, noise.sd = .1
  )
    
  #---- Apply methods ----
  results <- list()
  results[["CART"]] <- apply(dataSim$Ysim, 2, CART.apply, xb = dataSim$Xsim)
  results[["MARS"]] <- apply(dataSim$Ysim, 2, MARS.apply, 
    xb = dataSim$Xsim, p = p)
  results[["PRIM"]] <- apply(dataSim$Ysim, 2, PRIM.apply, xb = dataSim$Xsim)
#  results[["Classical"]] <- apply(dataSim$Ysim, 2, classical.apply, 
#    xb = dataSim$Xsim, extremes = dataSim$extPred > 0)
  
  #---- Result comparison ----
  
  # All thresholds estimated
  thresholds <- lapply(results, sapply, "[[", "thresholds")
  sTrue <- mapply("quantile", as.data.frame(dataSim$Xsim), s)
    
  # Bias
  meanThresh <- sapply(thresholds, apply, 1, mean, na.rm = T)
  bias <- meanThresh - matrix(sTrue, nrow = p, ncol = length(results))
  biasRes[[i]] <- bias
  
  # Standard error
  SE <- sapply(thresholds, apply, 1, sd, na.rm = T)
  SEres[[i]] <- SE
  
  # Mean Square Error
  sDiff <- sapply(thresholds, "-", matrix(sTrue, nrow = p, ncol = B), 
    simplify = "array")
  RMSE <- apply(sDiff, c(1,3), function(x) sqrt(sum(x^2, na.rm = T)))
  RMSEres[[i]] <- RMSE
  
  #---- Saving Results ----
  if (SAVE){
    name <- do.call(paste, c(sapply(combParam[i,], names), list(sep = "_")))
    out.dir <- sprintf("%s/%s", savepath, name)
    if(!dir.exists(out.dir)) dir.create(out.dir, recursive = T)
    setwd(out.dir)
    
    # Results
    export_list(dataSim$parameters, file = "Parameters.txt")
    save(results, bias, RMSE, file = "All_results.RData")
    
    # Plots
    # Plotting the true relationships
    if (L > 0){
      png(filename = "True_relationships.png")
      par(mfrow = n2mfrow(p), mar = c(1, 1, 3, 1))
      for (j in 1:p){
        persp(x = seq(0, 1, length.out = 20), y = 0:L, 
          dataSim$fTrue[[j]], xlab = "x", ylab = "Lag", zlab = "f", 
          main = colnames(dataSim$Xsim)[j],
          col = "lightskyblue", theta = 230, phi = 30, ticktype = "detailed")
      }
      dev.off()
    }
    
    Xints <- apply(dataSim$Xsim, 2, cut, 20, labels = FALSE)
    surf_cut <- aggregate(dataSim$linPred + dataSim$extPred, 
      by = as.data.frame(Xints), mean, na.rm = T)
    colnames(surf_cut)[3] <- "Y"
    
    png(filename = "TrueBivariate.png")
    wireframe(Y ~ V1 + V2, data = surf_cut, drape = T)
    dev.off()
    
    # Response illustration
    png(filename = "Realization1_dim2.png")
    plot(dataSim$Xsim, pch = 16, 
      col = heat.colors(10)[cut(dataSim$Ysim[,1], 10)], 
      xlab = "X1", ylab = "X2")
    rect(s, s, par("usr")[2] + 1,  par("usr")[4] + 1, 
      border = "black", lwd = 3, lty = 2)
    dev.off()
      
    png(filename = "Realization1_dim1.png")
    par(mfrow = n2mfrow(p))
    for (j in 1:p){
      plot(dataSim$Xsim[,j], dataSim$Ysim[,1], pch = 16, 
        xlab = sprintf("X%i", j), ylab = "Y")
      abline(v = s, lwd = 3, lty = 2, col = "red")
    }
    dev.off()
    
    png(filename = "Realization1_ACF.png")
    par(mfrow = c(2,1))
    acf(dataSim$Ysim[,1], na.action = na.pass, main = "")
    pacf(dataSim$Ysim[,1], na.action = na.pass, main = "")
    dev.off()
    
    # Thresholds
    png(filename = "Thresholds.png")
    par(mfrow = n2mfrow(p))
    for (j in 1:p){
      jthresh <- sapply(thresholds, "[", j, )
      boxplot(jthresh, border = cols, lwd = 2, ylab = "Threshold", xlab = "",
        main = colnames(dataSim$Xsim)[j], cex.lab = 1.3, cex.axis = 1.2, 
        varwidth = T)
      abline(h = s, lwd = 3, lty = 2)
    }
    dev.off()
    
    # Criteria
    png(filename = "Bias.png")
    bp <- barplot(t(bias), col = cols, border = NA, ylab = "Bias", 
      beside = TRUE, 
      main = colnames(dataSim$Xsim)[j], cex.lab = 1.3, cex.axis = 1.2, 
      cex.names = 1.2, 
      ylim = range(c(bias, 0)) + diff(range(c(bias, 0))) * c(-.1, .1))
    abline(h = 0)
    text(bp, t(bias), formatC(t(bias), format = "e", digits = 0), 
      pos = 2 + sign(t(bias)), cex = 1.2, xpd = T)
    outerLegend("topcenter", names(thresholds), fill = cols, bty = "n", 
      ncol = 3)
    dev.off()
    
    png(filename = "RMSE.png")
    bp <- barplot(t(RMSE), col = cols, border = NA, ylab = "RMSE", 
      beside = TRUE, 
      main = colnames(dataSim$Xsim)[j], cex.lab = 1.3, cex.axis = 1.2, 
      cex.names = 1.2,
      ylim = range(c(RMSE, 0)) + diff(range(c(RMSE, 0))) * c(-.1, .1))
    abline(h = 0)
    text(bp, t(RMSE), formatC(t(RMSE), format = "e", digits = 0), 
      pos = 2 + sign(t(RMSE)), cex = 1.2, xpd = T)
    outerLegend("topcenter", names(thresholds), fill = cols, bty = "n", 
      ncol = 3)
    dev.off()
  }
}

if (SAVE){ 
  save(combParam, biasRes, SEres, RMSEres, 
    sprintf("%s/AllResults.RData", savepath))
}

##############################################################################
#
#                    Plotting the results

source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Internet/Nested_loop_plot.r")

load(sprintf("%s/AllResults.RData", savepath))

respath <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/ArticleV5"

# Parameters
varnames <- c("extType", "extBetas", "Lsim")
varlabs <- c("fs", expression(beta), "L")

# Prepare data generating parameters
combs <- as.data.frame(lapply(combParam, names))
combs <- within(combs, {
  extType <- as.integer(extType)
  extBetas <- as.integer(substring(extBetas, 2)) / 100
  Lsim <- as.integer(substring(Lsim, 2))
#  Xdep <- as.numeric(substring(Xdep, 5)) / 10
})

#---------------------------------------------------
#             Nested loop plots
#---------------------------------------------------

### Bias
png(filename = sprintf("%s/Simulations_NestedPlot_Bias.png", respath), 
  units = "in", res = 100, height = 10, width = 7)

par(mar = c(2, 4, 3, 2) + .1)
layoutmat <- rbind(5, 1:2, 6, 3:4, 7)
layout(layoutmat, height = c(.1, .9, .1, .9, .1), width = c(.8, .2))

# Loop on plots
for (j in 1:p){
  biasj <- sapply(biasRes, "[", j, )
  resj <- data.frame(combs, t(biasj))
  
  pllim <- range(biasj)
  lilim <- pllim[1] - c(.5, .1) * diff(pllim)
    
  # Reorder for Nested loop plot
  nldata <- nestedloop(resj, varnames = varnames,
    varlabels= varlabs)
  nldata$extType <- factor(nldata$extType, levels = c(1, 2), 
    labels = c("Jump", "Break"))

  # Finally plot    
  matplot(nldata[rownames(biasj)], type = "n", cex.lab = 1.3,   
    xlab = "", ylab = "Bias", las = 1, xaxt = "n", 
    ylim = c(lilim[1], pllim[2]))
  lines(nldata, lty = 3, col = "lightgrey")
  lines(nldata, which = "r", ymin.refline = lilim[1], ymax.refline = lilim[2], 
    log = F, col = "black")
  matlines(nldata[rownames(biasj)], type = "s", col = cols, lty = 1)
  abline(h = 0, lwd = 2)
  title(main = "Detail")
#    outerLegend("topcenter", rownames(biasj), lwd = 2, col = cols, ncol = 3,
#      bty = "n")
  
  overallBias <- apply(resj[,rownames(biasj)], 2, mean)
#  par(mar = c(2, 0, 3, 1) + .1)
  barplot(overallBias, col = cols, xlab = "", ylab = "", 
    names.arg = "", main = "Overall")
}
par(mar = rep(0,4))
plot.new()
text(x = .5, y = .5, labels = "a) X1", cex = 2)
plot.new()
text(x = .5, y = .5, labels = "b) X2", cex = 2)
plot.new()
legend("center", rownames(biasj), lwd = 2, col = cols, ncol = 3, bty = "n",
  cex = 1.5)
  
dev.off()
  
### Standard error
png(filename = sprintf("%s/Simulations_NestedPlot_SE.png", respath), 
  units = "in", res = 100, height = 10, width = 7)

par(mar = c(2, 4, 3, 2) + .1)
layoutmat <- rbind(5, 1:2, 6, 3:4, 7)
layout(layoutmat, height = c(.1, .9, .1, .9, .1), width = c(.8, .2))

# Loop on plots
for (j in 1:p){
  SEj <- sapply(SEres, "[", j, )
  resj <- data.frame(combs, t(SEj))
  
  pllim <- range(SEj)
  lilim <- pllim[1] - c(.5, .1) * diff(pllim)
    
  # Reorder for Nested loop plot
  nldata <- nestedloop(resj, varnames = varnames,
    varlabels= varlabs)
  nldata$extType <- factor(nldata$extType, levels = c(1, 2), 
    labels = c("Jump", "Break"))

  # Finally plot    
  matplot(nldata[rownames(SEj)], type = "n", cex.lab = 1.3,   
    xlab = "", ylab = "Standard error", las = 1, xaxt = "n", 
    ylim = c(lilim[1], pllim[2]))
  lines(nldata, lty = 3, col = "lightgrey")
  lines(nldata, which = "r", ymin.refline = lilim[1], ymax.refline = lilim[2], 
    log = F, col = "black")
  matlines(nldata[rownames(SEj)], type = "s", col = cols, lty = 1)
  abline(h = 0, lwd = 2)
  title(main = "Detail")
#    outerLegend("topcenter", rownames(biasj), lwd = 2, col = cols, ncol = 3,
#      bty = "n")
  
  overallSE <- apply(resj[,rownames(SEj)], 2, mean)
#  par(mar = c(2, 0, 3, 1) + .1)
  barplot(overallSE, col = cols, xlab = "", ylab = "", 
    names.arg = "", main = "Overall")
}
par(mar = rep(0,4))
plot.new()
text(x = .5, y = .5, labels = "a) X1", cex = 2)
plot.new()
text(x = .5, y = .5, labels = "b) X2", cex = 2)
plot.new()
legend("center", rownames(SEj), lwd = 2, col = cols, ncol = 3, bty = "n",
  cex = 1.5)
  
dev.off()

  

##############################################################################
#
#                    Illustrative plots

# generate data to obtain the true lag-response functions
gener <- generate.data(B = B, Lsim = 5, p = 2, ffuns = "flinear", 
  fbetas = list(c(0,1), c(0,1)), wfuns = c("wconstant", "wdecay"), 
  wbetas = list(1 / 6, c(1, -5 / 6))
)

lagTrue <- lapply(gener$fTrue, "[", 20, )
lagTrue <- lapply(lagTrue, c, rep(0,5))

# Generate other data to obtain true impact of thresholds
generJump <- generate.data(B = 1, Lsim = 1, p = 1, ffuns = "flinear", 
  fbetas = list(c(0,1)), wfuns = "wone", wbetas = list(1), s = 0.8,
  extBetas = c(.5, 0) 
)
fJump <- as.data.frame(generJump[c("Xsim", "linPred", "extPred")])
fJump <- within(fJump, f <- linPred + extPred)
fJump <- fJump[order(fJump$Xsim),]



generBreak <- generate.data(B = 1, Lsim = 1, p = 1, ffuns = "flinear", 
  fbetas = list(c(0,1)), wfuns = "wone", wbetas = list(1), s = 0.8,
  extBetas = c(0, 2) 
)
fBreak <- as.data.frame(generBreak[c("Xsim", "linPred", "extPred")])
fBreak <- within(fBreak, f <- linPred + extPred)
fBreak <- fBreak[order(fBreak$Xsim),]

# Plot
png(filename = "Simulation_designs.png", width = 7, height = 10, res = 100, 
  units = "in")

par(mfrow = c(2,1), mar = c(5, 5, 4, 2) + .1)
plot(0:10, lagTrue[[2]], xlab = "Lag", ylab = "f", type = "b", lwd = 2,
  pch = 16, cex.axis = 1.2, cex.lab = 1.3, col = "darkgrey", 
  main = "a) Lag functions")
lines(0:10, lagTrue[[1]], type = "b", pch = 17, lty = 2, col = "black", lwd = 2)
legend("topright", c(expression(x[1]), expression(x[2])), 
  col = c("black", "darkgrey"), lwd = 2, lty = 2:1, pch = 17:16, bty = "n")

plot(fJump[,c("Xsim", "f")], xlab = "x", ylab = "f", type = "l", lwd = 2,
  cex.axis = 1.2, cex.lab = 1.3, main = "b) Extreme designs")
lines(fBreak[,c("Xsim", "f")], lwd = 2, col = "darkgrey", lty = 2)
legend("topleft", c("Jump", "Break"), col = c("black", "darkgrey"), lty = 1:2,
  lwd = 2, bty = "n")

dev.off()