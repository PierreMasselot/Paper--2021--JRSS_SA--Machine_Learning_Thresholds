#########################################################
#
#            Parametric simulations in paper
#                                                  
#########################################################

library(quantmod)
library(earth)
library(hhws)
library(lattice)
library(primr)
library(mgcv)
library(partykit)
library(AIM)
library(parallel)

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/1 - Thresholds/Paper--Machine-Learning-Thresholds")

source("0 - Misc_functions.R")
source("0 - Threshold_functions.R")

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

SAVE <- FALSE
respath <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/Article V8"

# Important constant parameters
p <- 2
s <- .95  # Quantile
B <- 10

# Varying parameters between simulations
varParams <- list()
#  varParams$Lsim <- 0:5; names(varParams$Lsim) <- sprintf("L%s", Lsim)
varParams$extType <- list(jump = c(1, 0, 0), linear = c(0, 1, 1))
varParams$extBetas <- list(X10 = .1, X30 = .3, X50 = .5, X100 = 1, 
                   X300 = 3, X500 = 5, X1000 = 10)


combParam <- do.call(expand.grid, varParams)
nc <- nrow(combParam)
nvp <- ncol(combParam)

# Other useful objects
cols <- c("forestgreen", "cornflowerblue", "firebrick", "goldenrod", 
  "slategrey")

# Scores of interest 
biasRes <- vector("list", nc)
SEres <- vector("list", nc)
RMSEres <- vector("list", nc)

# Uncertainty
biasSE <- SESE <- RMSESE <- vector("list", nc)

#---- Paralell computing
# Initialize cluster for parallel computation
cl <- makeCluster(2)
# Transfer objects in cluster
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(quantmod)
  library(earth)
  library(hhws)
  library(lattice)
  library(primr)
  library(mgcv)
  library(partykit)
  library(AIM)
})

#---- Loop for simulations ---
for(i in 1:nc){  
  print(sprintf("%i / %i", i, nc)); flush.console()
  
  #---- Generate data
  L <- 3
  XBetas <- combParam[[i, "extType"]] * combParam[[i, "extBetas"]]
  
  dataSim <- generate.data(B = B, Lsim = L, p = p, s = s,
    obetas = c(0, 1, 1), ffuns = "flinear", fbetas = list(c(0,1), c(0,1)), 
    wfuns = c("wconstant", "wdecay"), 
    wbetas = list(1 / (L + 1), c(1, -5 / (L + 1))),
    extBetas = XBetas, XdepOrder = .6, noise.sd = .1
  )
  
  # Transfer objects in cluster
  clusterExport(cl, "dataSim")
    
  #---- Apply methods ----
  results <- list()
  results[["MOB"]] <- parApply(cl, dataSim$Ysim, 2, MOB.apply, xb = dataSim$Xsim)
  results[["MARS"]] <- parApply(cl, dataSim$Ysim, 2, MARS.apply, 
    xb = dataSim$Xsim, degree = p)
  results[["PRIM"]] <- parApply(cl, dataSim$Ysim, 2, PRIM.apply, xb = dataSim$Xsim)
  results[["AIM"]] <- parApply(cl, dataSim$Ysim, 2, AIM.apply, xb = dataSim$Xsim, 
    backfit = T, numcut = 1)
  results[["GAM"]] <- parApply(cl, dataSim$Ysim, 2, GAM.apply, xb = dataSim$Xsim)
  
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
}

if (SAVE){ 
  save(combParam, biasRes, SEres, RMSEres, 
    file = sprintf("%s/Results_simulations.RData", respath))
}

##############################################################################
#
#                    Plotting the results

load(sprintf("%s/Results_simulations.RData", respath))

# Parameters
varnames <- c("extType", "extBetas")
varlabs <- c("fs", expression(beta))

# Prepare data generating parameters
combs <- as.data.frame(lapply(combParam, names))
combs <- within(combs, {
  extType <- factor(as.integer(extType), levels = c(1, 2), 
    labels = c("Jump", "Break"))
  extBetas <- as.integer(substring(extBetas, 2)) / 100
})

# Number of models
nm <- ncol(RMSEres[[1]])

#----------------------------------------------------------------
#                         Mean plot
#----------------------------------------------------------------
  
ylim <- range(unlist(RMSEres))
  
x11(width = 15)
layout(rbind(matrix(1:4, nrow = 2, byrow = T), 5), height = c(0.45, 0.45, 0.1))
par(mar = c(2, 2, 4, 2) + .1, oma = c(4, 5, 2, 0),
  cex.lab = 1.5, cex.main = 2, cex.axis = 1.3, xpd = NA)
for (j in 1:p){
  # Extract and format data
  jdat <- as.data.frame(t(sapply(RMSEres, "[", j, )))
  jdat <- as.data.frame(split(jdat, combs$extType))
  xsim <- unique(combs[,"extBetas"])
  
  # First extreme type
  matplot(jdat[,1:nm], xaxt = "n", ylim = ylim, yaxt = "n", 
    type = "b", col = cols, lwd = 2, pch = 14 + 1:nm, 
    cex = 1.5, log = "y", ylab = "RMSE (log)",
    main = bquote(.(letters[j * 2 - 1]) * ")" ~ X[.(j)] * "," ~ 
      .(as.character(unique(combs$extType)[1]))))
  axis(1, at = seq_along(xsim), labels = xsim, cex.axis = 1.2)
  if (j == p) title(xlab = varlabs[2], xpd = NA)
  axis(2, at = axTicks(2, axp = c(.1, 3, 1), log = T))
  
  # Second extreme type
  matplot(jdat[,nm + 1:nm], xaxt = "n", ylim = ylim, 
    type = "b", col = cols, lwd = 2, pch = 14 + 1:nm, 
    ylab = "", yaxt = "n", 
    cex = 1.5, log = "y",
    main = bquote(.(letters[j * 2]) * ")" ~ X[.(j)] * "," ~ 
      .(as.character(unique(combs$extType)[2]))))
  axis(1, at = seq_along(xsim), labels = xsim, cex.axis = 1.2)
  if (j == p) title(xlab = varlabs[2], cex.lab = 1.5, xpd = NA)
  axis(2, at = axTicks(2, axp = c(.1, 3, 1), log = T))    
}
# Legend in a new panel
par(mar = rep(0,4))
plot.new()
legend("bottom", colnames(RMSEres[[1]]), col = cols, pch = 14 + 1:nm, 
  lty = 1:5, ncol = nm, cex = 1.5)

if (SAVE){
  dev.print(png, filename = sprintf("%s/Fig2_SimulationResults.png", respath), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/Fig2_SimulationResults.eps", respath))
}


#----------------------------------------------------------------
#       Illustration of the simulation designs
#----------------------------------------------------------------

# generate data to obtain the true lag-response functions
gener <- generate.data(B = 1, Lsim = 3, p = 2, ffuns = "flinear", 
  fbetas = list(c(0,1), c(0,1)), wfuns = c("wconstant", "wdecay"), 
  wbetas = list(1 / 6, c(1, -5 / 6))
)

lagTrue <- lapply(gener$fTrue, "[", 20, )
lagTrue <- lapply(lagTrue, c, rep(0,2))

# Generate other data to obtain true impact of thresholds
generJump <- generate.data(B = 1, Lsim = 1, p = 1, ffuns = "flinear", 
  fbetas = list(c(0,1)), wfuns = "wone", wbetas = list(1), s = 0.8,
  extBetas = c(1, 0) 
)
fJump <- as.data.frame(generJump[c("Xsim", "linPred", "extPred")])
fJump <- within(fJump, f <- linPred + extPred)
fJump <- fJump[order(fJump$Xsim),]



generBreak <- generate.data(B = 1, Lsim = 1, p = 1, ffuns = "flinear", 
  fbetas = list(c(0,1)), wfuns = "wone", wbetas = list(1), s = 0.8,
  extBetas = c(0, 5) 
)
fBreak <- as.data.frame(generBreak[c("Xsim", "linPred", "extPred")])
fBreak <- within(fBreak, f <- linPred + extPred)
fBreak <- fBreak[order(fBreak$Xsim),]

# Plot
x11(height = 10)
par(mfrow = c(2,1), mar = c(5, 5, 4, 2) + .1)
plot(0:5, lagTrue[[2]], xlab = "Lag", ylab = "f", type = "b", lwd = 2,
  pch = 16, cex.axis = 1.2, cex.lab = 1.3, col = "darkgrey", 
  main = "a) Lag functions")
lines(0:5, lagTrue[[1]], type = "b", pch = 17, lty = 2, col = "black", lwd = 2)
legend("topright", c(expression(x[1]), expression(x[2])), 
  col = c("black", "darkgrey"), lwd = 2, lty = 2:1, pch = 17:16, bty = "n")

plot(fJump[,c("Xsim", "f")], xlab = "x", ylab = "f", type = "l", lwd = 2,
  cex.axis = 1.2, cex.lab = 1.3, main = "b) Extreme designs")
lines(fBreak[,c("Xsim", "f")], lwd = 2, col = "darkgrey", lty = 2)
legend("topleft", c("Jump", "Break"), col = c("black", "darkgrey"), lty = 1:2,
  lwd = 2, bty = "n")

if (SAVE){
  dev.print(png, filename = sprintf("%s/Fig1_SimulationDesign.png", respath), 
    units = "in", res = 100)
  dev.copy2eps(file = sprintf("%s/Fig1_SimulationDesign.eps", respath))
}