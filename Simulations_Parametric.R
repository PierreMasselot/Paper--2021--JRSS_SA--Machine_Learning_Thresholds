#########################################################
#
#            Parametric simulations in paper
#                                                  
#########################################################

library(quantmod)
library(rpart)
library(earth)
library(hhws)

source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/PRIM/R/PRIM_functions.R")
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
#' @param extBetas Parameters for the linear function of extremes.
#' @param XdepOrder AR order for temporal dependence of indicators.
#' @param YdepOrder AR ordre for temporal dependence of response.
#' @param noise.sd Standard deviation of the random part of response, relative
#'    to the standard deviation of the deterministic part.
generate.data <- function(B = 1000, n = 5000, p = 2, Lsim = 1, 
  obetas = 1, ffuns = "fconstant", fbetas = list(), wfuns = "wone", 
  wbetas = list(), s = .8, extBetas = 1, XdepOrder = 0, YdepOrder = 0,
  noise.sd = .2)
{
  params <- as.list(environment())
  
  # Recycling necessary parameters
  ffuns <- rep_len(ffuns, p)
  wfuns <- rep_len(wfuns, p)
  s <- rep_len(s, p)
  obetas <- rep_len(obetas, p + 1)
  wbetas <- rep_len(wbetas, p + 1)
  
  # Predictors X
  Xsim <- replicate(p, arima.sim(list(ar = XdepOrder), n, runif))
  Xsim <- apply(Xsim, 2, fscaling)
  
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
  
  # Extreme part
  uni.extremes <- mapply(">=", as.data.frame(Xsim), s)
  extremes <- apply(uni.extremes, 1, all)
  
  extBetas <- extBetas * diff(range(linPred, na.rm = T))
  Xext <- mapply("-", as.data.frame(Xsim), s)
  extPred <- ifelse(extremes, cbind(1, Xext) %*% extBetas, 0)

  # Response simulation
  Ysim <- replicate(B, linPred + extPred + 
    arima.sim(list(ar = YdepOrder), n, rnorm, sd = noise.sd), simplify = TRUE)
  
  output <- list(Ysim = Ysim, Xsim = Xsim, linPred = linPred, extPred = extPred,
    fTrue = fTrue, overall = overall, parameters = params)
  return(output)
}


#---------------------------------------------------
#             Threshold finding methods
#---------------------------------------------------

#' CART 
CART.apply <- function(yb, xb){
  datab <- data.frame(Y = yb, xb)
  # Grow the tree
  treeb <- rpart(Y ~ ., data = datab, method = "anova", 
    control = rpart.control(minsplit = 10, cp = 0.0001))
  minind <- which.min(treeb$cptable[,"xerror"])
  best_cp <- treeb$cptable[minind, "CP"]
  treeb <- prune(treeb, best_cp)
  
  # Extract thresholds
  boxb <- get.box(treeb)
  
  # Alarms
  uni.alb <- mapply(">=", datab[,-1], boxb$box[1,])
  alarmb <- apply(uni.alb, 1, all)
   
  list(thresholds = boxb$box[1,], alarms = datab[alarmb,1])
}

#' MARS
#'
#' @param p Interaction degree
MARS.apply <- function(yb, xb, p = 2){
  datab <- data.frame(Y = yb, xb)
  # Apply MARS
  marsb <- earth(Y ~ ., data = na.omit(datab), degree = 2)
  
  # Extract alarms
  uni.alb <- mapply(">=", datab[,-1], 
    apply(marsb$cuts, 2, max))
  alarmb <- apply(uni.alb, 1, all) 
  
  list(thresholds = apply(marsb$cuts, 2, max), alarms = datab[alarmb,1])
}

#' PRIM
PRIM.apply <- function(yb, xb){
  xb <- xb[complete.cases(yb),]
  yb <- na.omit(yb)
  n <- length(yb)
  
  # Peeling the box
  peelb <- peeling.sequence(yb, xb, alpha = .05, beta.stop = 6/n, 
    peeling.side = "left")
  
  # Extract box : box following the biggest increase in mean
  supb <- peelb$support[which.max(diff(peelb$yfun))]
  chosenb <- extract.box(peelb, supb)
  boxb <- pasting.sequence(yb, xb, 
    small.box = chosenb$limits, peeling.side = "left", alpha = 0.01)
    
  list(thresholds = boxb$limits[1,], 
    alarms = in.box(xb, boxb$limits, yb))
}

#' Classical
classical.apply <- function(yb, xb, extremes){
  xb <- xb[complete.cases(yb),]
  yb <- na.omit(yb)
  n <- length(yb)
  
  # Extract episodes
  OMt <- floor(quantile(yb, (n - sum(extremes)) / n) / 5) * 5
  OM.episodes <- episodes(yb, OMt, r = 3, l = 3)
  
  # Tested thresholds
  tested <- find.threshold(as.data.frame(xb), episodes = OM.episodes, 
    u.grid = seq(.5, 1, by = .05), thinning = "none")
  youden <- tested[,"Sensitivity"] + tested[,"Specificity"] - 1
  
  # Thresholds
  thresholds <- tested[which.max(youden),grep("threshold", names(tested))]
    
  # Extract alarms
  uni.alb <- mapply(">=", as.data.frame(xb), thresholds)
  alarmb <- apply(uni.alb, 1, all)
  
  # Result
  list(thresholds = unlist(thresholds), alarms = yb[alarmb]) 
}

##############################################################################
#
#                            Simulations 

sim.name <- "Test"
SAVE <- FALSE

# Important constant parameters
p <- 2
s <- .8
B <- 2

# Varying parameters between simulations
Lsim <- 1:5; names(Lsim) <- sprintf("L%s", Lsim)
extBetas <- list(X10 = c(0,.1,.1), X30 = c(0,.3,.3), X50 = c(0,.5,.5), 
  X100 = c(0,1,1), X300 = c(0,3,3), X500 = c(0,5,5), X1000 = c(0,10,10))

combParam <- expand.grid(1:length(Lsim), 1:length(extBetas))
nc <- nrow(combParam)

# Other useful object
cols <- c("forestgreen", "cornflowerblue", "skyblue", "firebrick", "darkgrey")

#---- Loop for simulations ---=
for(i in 1:nc){  
  print(sprintf("%i / %i", i, nc)); flush.console()
  
  #---- Generate data
  dataSim <- generate.data(B = B, Lsim = Lsim[combParam[i,1]], p = p, s = s,
    obetas = c(0, 1, 1), ffuns = "flinear", fbetas = list(c(0,1), c(0,1)), 
    wfuns = c("wconstant", "wdecay"), 
    wbetas = list(1 / (Lsim[combParam[i,1]] + 1), 
      c(1, -5 / Lsim[combParam[i,1]])),
    extBetas = extBetas[[combParam[i,2]]], XdepOrder = .6, noise.sd = .1
  )
    
  #---- Apply methods ----
  results <- list()
  results[["CART"]] <- apply(dataSim$Ysim, 2, CART.apply, xb = dataSim$Xsim)
  results[["MARS"]] <- apply(dataSim$Ysim, 2, MARS.apply, xb = dataSim$Xsim)
  results[["MARS_inter"]] <- apply(dataSim$Ysim, 2, MARS.apply, 
    xb = dataSim$Xsim, p = p)
  results[["PRIM"]] <- apply(dataSim$Ysim, 2, PRIM.apply, xb = dataSim$Xsim)
  results[["Classical"]] <- apply(dataSim$Ysim, 2, classical.apply, 
    xb = dataSim$Xsim, extremes = dataSim$extPred > 0)
  
  #---- Result comparison ----
  
  # All thresholds estimated
  thresholds <- lapply(results, sapply, "[[", "thresholds")
    
  # Bias
  meanThresh <- sapply(thresholds, apply, 1, mean, na.rm = T)
  bias <- meanThresh - matrix(s, nrow = p, ncol = length(results))
  
  # Mean Square Error
  sDiff <- sapply(thresholds, "-", matrix(s, nrow = p, ncol = B), 
    simplify = "array")
  RMSE <- apply(sDiff, c(1,3), function(x) sqrt(sum(x^2, na.rm = T)))
  
  #---- Saving Results ----
  if (SAVE){
    name <- paste(sim.name, names(Lsim)[combParam[i,1]], 
      names(extBetas)[combParam[i,2]], sep = "_")
    out.dir <- sprintf("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/Simulations/%s", sim.name)
    if(!dir.exists(out.dir)) dir.create(out.dir, recursive = T)
    setwd(out.dir)
    
    # Results
    export_list(dataSim$parameters, file = "Parameters.txt")
    save(results, bias, RMSE, file = "All_results.RData")
    
    # Plots
    # Plotting the true relationships
    png(filename = "True_relationships.png")
    par(mfrow = n2mfrow(p), mar = c(1, 1, 3, 1))
    for (j in 1:p){
      persp(x = seq(0, 1, length.out = 20), y = 0:combParam[[1]][i], 
        dataSim$fTrue[[j]], xlab = "x", ylab = "Lag", zlab = "f", 
        main = colnames(dataSim$Xsim)[j],
        col = "lightskyblue", theta = 230, phi = 30, ticktype = "detailed")
    }
    dev.off()
    
    Xints <- apply(dataSim$Xsim, 2, cut, 20, labels = FALSE)
    surf_cut <- aggregate(dataSim$linPred + dataSim$extPred, 
      by = as.data.frame(Xints), mean, na.rm = T)
    colnames(surf_cut)[3] <- "Y"
    
    png(filename = "TrueBivariate.png")
    lattice::wireframe(Y ~ V1 + V2, data = surf_cut, drape = T)
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

