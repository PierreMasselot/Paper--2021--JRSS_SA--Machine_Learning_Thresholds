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
#                 Initialization
#---------------------------------------------------

sim.name <- "JshapeLinear_X30_L0_noDep_s8"
SAVE <- FALSE

parameters <- within(list(),{
  # Parameters
  B <- 1000 # Number of simulations
  n <- 5000
  p <- 2
  Lsim <- 1
  
  # Function parameters
  obetas <- c(0, 1, 1) # Intercept is the first element
  ffuns <- c("fJshape", "flinear")
  fbetas <- list(c(0, 0.05, 0.5, 0.2, 0.1), c(0, .5))   
  wfuns <- c("wone", "wone")
  wbetas <- list(c(1, -5/Lsim), 1 / (Lsim + 1))
  
  s <- c(.8, .8)
  extBetas <- c(0, 10, 10)  # Coefficients for a linear function above the threshold.
                 # Relative to the range of Y
  
  XdepOrder <- 0
  YdepOrder <- 0
  
  noise.sd <- .2
})
attach(parameters)

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


# Object containing results
results <- list()
cols <- c("forestgreen", "cornflowerblue", "skyblue", "firebrick", "darkgrey")

#---------------------------------------------------
#                 Simulation
#---------------------------------------------------

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

# Plotting the true relationships
x11(title = "True_relationships")
par(mfrow = n2mfrow(p), mar = c(1, 1, 3, 1))
for (j in 1:p){
  persp(x = seq(0, 1, length.out = 20), y = 0:Lsim, fTrue[[j]], xlab = "x",
    ylab = "Lag", zlab = "f", main = colnames(Xsim)[j],
    col = "lightskyblue", theta = 230, phi = 30, ticktype = "detailed")
}

Xints <- apply(Xsim, 2, cut, 20, labels = FALSE)
surf_cut <- aggregate(linPred + extPred, by = as.data.frame(Xints), mean, 
  na.rm = T)
colnames(surf_cut)[3] <- "Y"

x11(title = "TrueBivariate")
lattice::wireframe(Y ~ V1 + V2, data = surf_cut, drape = T)

# Response simulation
Ysim <- replicate(B, linPred + extPred + 
  arima.sim(list(ar = YdepOrder), n, rnorm, sd = noise.sd), simplify = TRUE)

# Illustration response
x11(title = "Realization1_dim2")
plot(Xsim, pch = 16, col = heat.colors(10)[cut(Ysim[,1], 10)], xlab = "X1",
  ylab = "X2")
rect(s[1], s[2], par("usr")[2] + 1,  par("usr")[4] + 1, border = "black", 
  lwd = 3, lty = 2)
  
x11(title = "Realization1_dim1")
par(mfrow = n2mfrow(p))
for (j in 1:p){
  plot(Xsim[,j], Ysim[,1], pch = 16, xlab = sprintf("X%i", j), ylab = "Y")
  abline(v = s[j], lwd = 3, lty = 2, col = "red")
}

x11(title = "Realization1_ACF")
par(mfrow = c(2,1))
acf(Ysim[,1], na.action = na.pass, main = "")
pacf(Ysim[,1], na.action = na.pass, main = "")

#---------------------------------------------------
#                      CART
#---------------------------------------------------

CART.apply <- function(yb){
  datab <- data.frame(Y = yb, Xsim)
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

results[["CART"]] <- apply(Ysim, 2, CART.apply)

#---------------------------------------------------
#                     MARS
#---------------------------------------------------

MARS.apply <- function(yb){
  datab <- data.frame(Y = yb, Xsim)
  # Apply MARS
  marsb <- earth(Y ~ ., data = na.omit(datab))
  
  # Extract alarms
  uni.alb <- mapply(">=", datab[,-1], 
    apply(marsb$cuts, 2, max))
  alarmb <- apply(uni.alb, 1, all) 
  
  list(thresholds = apply(marsb$cuts, 2, max), alarms = datab[alarmb,1])
}

results[["MARS"]] <- apply(Ysim, 2, MARS.apply)


#---------------------------------------------------
#              MARS with interactions
#---------------------------------------------------

MARS.apply <- function(yb){
  datab <- data.frame(Y = yb, Xsim)
  # Apply MARS
  marsb <- earth(Y ~ ., data = na.omit(datab), degree = 2)
  
  # Extract alarms
  uni.alb <- mapply(">=", datab[,-1], 
    apply(marsb$cuts, 2, max))
  alarmb <- apply(uni.alb, 1, all) 
  
  list(thresholds = apply(marsb$cuts, 2, max), alarms = datab[alarmb,1])
}

results[["MARS_inter"]] <- apply(Ysim, 2, MARS.apply)


#---------------------------------------------------
#                     PRIM
#---------------------------------------------------

PRIM.apply <- function(yb){
  xb <- Xsim[complete.cases(yb),]
  yb <- na.omit(yb)
  
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

results[["PRIM"]] <- apply(Ysim, 2, PRIM.apply)


#---------------------------------------------------
#                     Classical
#---------------------------------------------------

classical.apply <- function(yb){
  xb <- Xsim[complete.cases(yb),]
  yb <- na.omit(yb)
  
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
  list(thresholds = thresholds, alarms = yb[alarmb]) 
}


result[["Classical"]] <- apply(Ysim, 2, classical.apply)

#---------------------------------------------------
#                 Results comparison
#---------------------------------------------------

thresholds <- lapply(results, sapply, "[[", "thresholds")

# All thresholds estimated
x11(title = "Thresholds")
par(mfrow = n2mfrow(p))
for (j in 1:p){
  jthresh <- sapply(thresholds, "[", j, )
  boxplot(jthresh, border = cols, lwd = 2, ylab = "Threshold", 
    main = colnames(Xsim)[j], cex.lab = 1.3, cex.axis = 1.2, xlab = "",
    varwidth = T)
  abline(h = s[j], lwd = 3, lty = 2)
}

# Bias
meanThresh <- sapply(thresholds, apply, 1, mean, na.rm = T)
bias <- meanThresh - matrix(s, nrow = p, ncol = length(results))

x11(title = "Bias")
bp <- barplot(t(bias), col = cols, border = NA, ylab = "Bias", beside = TRUE, 
  main = colnames(Xsim)[j], cex.lab = 1.3, cex.axis = 1.2, cex.names = 1.2,
  ylim = range(c(bias, 0)) + diff(range(c(bias, 0))) * c(-.1, .1))
abline(h = 0)
text(bp, t(bias), formatC(t(bias), format = "e", digits = 0), 
  pos = 2 + sign(t(bias)), cex = 1.2, xpd = T)
outerLegend("topcenter", names(thresholds), fill = cols, bty = "n", ncol = 3)


# Mean Square Error
sDiff <- sapply(thresholds, "-", matrix(s, nrow = p, ncol = B), 
  simplify = "array")
RMSE <- apply(sDiff, c(1,3), function(x) sqrt(sum(x^2, na.rm = T)))

x11(title = "RMSE")
bp <- barplot(t(RMSE), col = cols, border = NA, ylab = "RMSE", beside = TRUE, 
  main = colnames(Xsim)[j], cex.lab = 1.3, cex.axis = 1.2, cex.names = 1.2,
  ylim = range(c(RMSE, 0)) + diff(range(c(RMSE, 0))) * c(-.1, .1))
abline(h = 0)
text(bp, t(RMSE), formatC(t(RMSE), format = "e", digits = 0), 
  pos = 2 + sign(t(RMSE)), cex = 1.2, xpd = T)
outerLegend("topcenter", names(thresholds), fill = cols, bty = "n", ncol = 3)

#---------------------------------------------------
#                 Saving results
#---------------------------------------------------

if(SAVE){
  out.dir <- sprintf("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/Simulations/%s", sim.name)
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive = T)
  setwd(out.dir)

  graph.names <- c("True_relationships", "True_Bivariate", "Realization1_dim2", 
    "Realization1_dim1", "realization1_ACF", "Thresholds", "Bias", "RMSE")
  for (d in 1:length(dev.list())){
    dev.set(dev.list()[d])
    dev.print(png, sprintf("%s.png", graph.names[d]), units = "in",
      width = dev.size()[1], height = dev.size()[1], res = 100)
  }

  export_list(parameters, file = "Parameters.txt")

  save(results, bias, RMSE file = "All_results.RData")
}