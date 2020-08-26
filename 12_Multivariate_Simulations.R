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
library(segmented)
library(dlnm)

source("00_Misc_functions.R")
source("01_Threshold_functions.R")

source("10_Simulation_functions.R")

#----------------------------------
#  Parameters
#----------------------------------

# Saving final figures
SAVE <- FALSE
respath <- "Figures"

# Important constant parameters
n <- 5000 # Sample size
ne <- 15  # Number of extreme Y values generated
B <- 5  # Replication number
L <- 3

# Varying parameters between simulations
varParams <- list()
#  varParams$Lsim <- 0:5; names(varParams$Lsim) <- sprintf("L%s", Lsim)
varParams$extType <- list(jump = c(1, 0), linear = c(0, 1))
varParams$p <- list(1, 2, 3, 5, 10, 20)

# Graphical parameters
cols <- c("forestgreen", "cornflowerblue", "firebrick", "goldenrod", 
  "slategrey")
  
#----------------------------------
#  Preparation of result storing objects
#----------------------------------

# Number of cases considered
combParam <- do.call(expand.grid, varParams)
nc <- nrow(combParam)
nvp <- ncol(combParam)

# Scores of threshold estimation
biasRes <- SEres <- RMSEres <- vector("list", nc)

# Scores of extreme days prediction
sensitivity <- precision <- F1 <- F2 <- vector("list", nc)

#----------------------------------
#  Preparing parallel computing (socket, on windows)
#----------------------------------

# Initialize cluster 
cl <- makeCluster(2)

# Load the necessary packages on each cluster
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
  library(segmented)
  library(dlnm)
})

#----------------------------------
#  Loop on all cases
#----------------------------------

for(i in 1:nc){  
  print(sprintf("%i / %i", i, nc)); flush.console()
  
  #---- Prepare parameters
  # Beta parameters for the extreme function
  ibetas <- combParam[[i, "extType"]] 
  ip <- combParam[[i, "p"]] 
  XBetas <- c(ibetas[1], rep_len(ibetas[2], ip))
  
  # Compute threshold according to number of extreme values desired
  s <- 1 - (ne / n)^(1/ip)
  
  #---- Generate data
  dataSim <- generate.data(B = B, Lsim = L, p = ip, s = s, 
    obetas = c(0, rep(1, ip)), 
    ffuns = "flinear", fbetas = replicate(ip, c(0,1), simplify = FALSE), 
    wfuns = "wconstant", wbetas = replicate(ip, 1 / (L + 1), simplify = FALSE),
    extBetas = XBetas, XdepOrder = .6, noise.sd = .1
  )
  
  #---- Transfer objects in cluster
  clusterExport(cl, "dataSim")
    
  #---- Apply methods 
  results <- list()
  results[["MOB"]] <- parApply(cl, dataSim$Ysim, 2, MOB.apply, 
    xb = dataSim$Xsim, zb = dataSim$Xsim)
  results[["MARS"]] <- parApply(cl, dataSim$Ysim, 2, MARS.apply, 
    xb = dataSim$Xsim, degree = 2)
  results[["PRIM"]] <- parApply(cl, dataSim$Ysim, 2, PRIM.apply, 
    xb = dataSim$Xsim, zb = dataSim$Xsim)
  results[["AIM"]] <- parApply(cl, dataSim$Ysim, 2, AIM.apply, 
    xb = dataSim$Xsim, backfit = T, numcut = 1)
  results[["GAM"]] <- parApply(cl, dataSim$Ysim, 2, GAM.apply, 
    xb = dataSim$Xsim)
  results[["DLNM"]] <- parApply(cl, dataSim$Ysim, 2, DLNM.apply, 
    xb = dataSim$Xsim, 
    crosspars = replicate(ip, list(lag = L + 1, 
      argvar = list(fun = "ps", df = 10), 
      arglag = list(fun = "strata", df = 1)), simplify = FALSE
    )
  )
  results[["SEG"]] <- parApply(cl, dataSim$Ysim, 2, seg.apply, 
    xb = dataSim$Xsim, zb = dataSim$Xsim)
  
  #---- Error of threshold estimation
  # True thresholds  
  sTrue <- mapply("quantile", as.data.frame(dataSim$Xsim), s)
    
  # Bias
  if (ip == 1) results <- lapply(results, matrix, nrow = 1)
  meanThresh <- sapply(results, apply, 1, mean, na.rm = T)
  bias <- meanThresh - matrix(sTrue, nrow = ip, ncol = length(results))
  biasRes[[i]] <- bias
  
  # Standard error
  SE <- sapply(results, apply, 1, sd, na.rm = T)
  SEres[[i]] <- SE
  
  # Mean Square Error
  sDiff <- sapply(results, "-", matrix(sTrue, nrow = ip, ncol = B), 
    simplify = "array")
  RMSE <- apply(sDiff, c(1,3), function(x) sqrt(sum(x^2, na.rm = T)))
  RMSEres[[i]] <- RMSE
  
  #---- Ability to predict extreme days
  # True extreme days
  trueX <- which(dataSim$extPred != 0)
  ntrue <- length(trueX)
  
  # Compute scores for each simulation of each model
  scores <- lapply(results, function(tr){
    preddays <- mapply(extract_alarms, thresholds = as.data.frame(tr), 
      y = as.data.frame(dataSim$Ysim), 
      MoreArgs = list(x = dataSim$Xsim), SIMPLIFY = FALSE)
    preddays <- lapply(preddays, function(x) as.numeric(names(x)))
    found <- lapply(preddays, "%in%", x = trueX)
    sens <- sapply(found, sum) / ntrue
    prec <- sapply(found, sum) / sapply(preddays, length)
    prec[prec == Inf] <- 0
    F1 <- sapply(preddays, Fscore, trueX, 1)
    F2 <- sapply(preddays, Fscore, trueX, 2) 
    list(sensitivity = sens, precision = prec, F1 = F1, F2 = F2)
  })
  
  # Store results
  sensitivity[[i]] <- sapply(scores, "[[", "sensitivity")
  precision[[i]] <- sapply(scores, "[[", "precision") 
  F1[[i]] <- sapply(scores, "[[", "F1") 
  F2[[i]] <- sapply(scores, "[[", "F2") 
}

stopCluster(cl)

if (SAVE){ 
  save(combParam, biasRes, SEres, RMSEres, 
    file = sprintf("%s/11_Results_simulations_bivariate.RData", respath))
}

#----------------------------------
#  Plots
#----------------------------------

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