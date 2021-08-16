################################################################################
#                              R code 
#                                                  
#     Machine learning approaches to identify thresholds in a heat-health 
#                      warning system context
#         Journal of the Royal Statistical Society - Series A
#                               2021
#
#                        Simulation study
#
#                   Code Author: Pierre Masselot
#
################################################################################

# This script implements the simulation study of section 4.

#----- Packages
library(earth) # MARS
library(primr) # PRIM (custom package, must be installed from github 
                # devtools::install_github("PierreMasselot/primr"))
library(mgcv) # GAM
library(partykit) # MOB
library(AIM) # AIM
library(parallel) # For parallelisation
library(segmented) # Segmented regression

#----- Load other functions
# Miscellaneous functions
source("00_Misc_functions.R")
# Wrappers to extract thresholds through the proposed methods
source("01_Threshold_functions.R")
# Functions simulating data
source("10_Simulation_functions.R")

#----------------------------------
#  Parameters
#----------------------------------

#----- Important constant parameters
pext <- .015 # Probability of extreme
B <- 1000  # Replication number
n <- 1000 # Sample size

#----- Varying parameters between simulations
varParams <- list()

# Type of extreme effect
varParams$extType <- list(One = c(0, 1), Two = c(0, 1, 1), 
  Three = c(0, 1, 1, 1), Four = c(0, 1, 1, 1, 1))

# Magnitude of extremes
varParams$extBetas <- list(X10 = .1, X20 = .2, X30 = .3, X50 = .5, X100 = 1)

# Correlation between variables
varParams$rho <- c(0, .5)
  
#----------------------------------
#  Preparation of result storing objects
#----------------------------------

# Number of scenarios
combParam <- do.call(expand.grid, varParams)
nc <- nrow(combParam)
nvp <- ncol(combParam)

# Objects storing scores of interest
sensitivity <- precision <- F1 <- F2 <- vector("list", nc)

#----------------------------------
#  Preparing parallel computing (socket, on windows)
#----------------------------------

# Initialize cluster 
cl <- makeCluster(detectCores() - 2)

# Load the necessary packages on each cluster
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(earth)
  library(primr)
  library(mgcv)
  library(partykit)
  library(AIM)
  library(parallel)
  library(segmented)
})

#----------------------------------
#  Loop on all cases
#----------------------------------

for(i in 1:nc){  
  print(sprintf("%i / %i", i, nc)); flush.console()
  
  #---- Generate data
  XBetas <- combParam[[i, "extType"]] * combParam[[i, "extBetas"]]
  pi <- length(combParam[[i, "extType"]]) - 1
  # Compute thresholds according to number of and corr between vars
  #     Avoids curse-of-dimensionality issues
  Smat <- matrix(combParam[[i, "rho"]], nrow = pi, ncol = pi)
  diag(Smat) <- 1
  s <- qmvnorm(pext, tail = "upper.tail", mean = rep(0, pi), 
    sigma = Smat)$quantile
  
  # Simulate
  set.seed(12345 + i)
  dataSim <- replicate(B, 
    generate.data(n = n, p = pi, s = s, stype = "absolute",
      obetas = c(3, rep(1, pi)), ffuns = "flinear", fbetas = list(c(0,1)), 
      extBetas = XBetas, rho = combParam[[i, "rho"]], rand.gen = "rpois"
    ), simplify = F
  )
  
  #---- Transfer objects in cluster
  clusterExport(cl, "dataSim")
    
  #---- Apply methods 
  results <- list()
  results[["MOB"]] <- parLapply(cl, dataSim, function(x){
    MOB.apply(yb = x$Ysim, xb = x$Xsim, zb = x$Xsim, minsize = 5,
      family = "poisson")
  })
  results[["MARS"]] <- parLapply(cl, dataSim, function(x){
    MARS.apply(yb = x$Ysim, xb = x$Xsim, degree = 2, endspan = 5,
      glm = list(family = "poisson"))
  })
  results[["PRIM"]] <- parLapply(cl, dataSim, function(x){
    PRIM.apply(yb = x$Ysim, xb = x$Xsim, zb = x$Xsim, beta.stop = 5/n,
      peeling.side = -1, family = "poisson")
  })
  results[["AIM"]] <- parLapply(cl, dataSim, function(x){
    AIM.apply(yb = x$Ysim, xb = x$Xsim, backfit = T, numcut = 3, 
      mincut = 5/n)
  })
  results[["SEG"]] <- parLapply(cl, dataSim, function(x) {
    seg.apply(yb = x$Ysim, xb = x$Xsim, zb = x$Xsim, 
      segpars = list(psi = list(NA), 
        control = seg.control(alpha = 5 / n, fix.npsi = F, n.boot = 0)),
      glmpars = list(family = "poisson"))
  })
  results[["GAM"]] <- parLapply(cl, dataSim, function(x){
    GAM.apply(yb = x$Ysim, xb = x$Xsim, family = "poisson")
  })

  #---- Ability to predict extreme days
  # True extreme days from simulated data
  trueX <- lapply(dataSim, function(x) which(x$extremes))
  ntrue <- sapply(trueX, length)

  # Compute scores for each simulation of each model
  scores <- lapply(results, function(x){
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
  sensitivity[[i]] <- sapply(scores, "[[", "sensitivity")
  precision[[i]] <- sapply(scores, "[[", "precision")
  F1[[i]] <- sapply(scores, "[[", "F1")
  F2[[i]] <- sapply(scores, "[[", "F2")
}

stopCluster(cl)

save.image(file = "Results/11_Results_simulations.RData")

