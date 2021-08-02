#########################################################
#
#            Parametric simulations
#         Part 0: Simulation functions
#                                                  
#########################################################

#---------------------------------------------------
#                 Relationship functions
#---------------------------------------------------

#----- Exposure-response functions

# Linear effect
#   This is the one used in the final simulation version
flinear <- function(x, Betas) Betas[1] + Betas[2] * x

# Constant (no effect of exposure) - Unused in the final manuscript
#   Unused in final manuscript
fconstant <- function(x, Betas) Betas

# J-shape mimicking temperature-mortality associations
#   Unused in final manuscript
fJshape <- function(x, Betas){
  pow <- outer(x, 0:4, "^")
  pow %*% Betas
}

#---------------------------------------------------
#              Data generating function
#---------------------------------------------------

# This function generates simulated data to be used for threshold evaluation

# Parameters
#   n   Number of simulated observations
#   p   Number of exposure variables
#   ncat  Number of categories for categorical variables. NA for continuous
#         variables
#   obetas  Outer betas: relative weights of the different exposure variables
#   ffuns   Character vector of exposure-response functions
#   fbetas  List of parameters for each function in ffuns
#   s   Thresholds
#   stype   How s is defined, either as absolute value, either as quantile of
#           of the distribution of its respective exposure variables
#   extBetas    Magnitude of extreme effect
#   YdepOrder   Autoregressive order for simulated response. 
#               0 for no autocorrelation
#   rho   Correlation between exposures
#   rand.gen    Random generator: either rnorm or rpois
#   noise.sd    Noise standard deviation

generate.data <- function(n = 1000, p = 2, ncat = NA, 
  obetas = 1, ffuns = "fconstant", fbetas = list(), 
  s = .8, stype = c("quantile", "absolute"), extBetas = 1,
  YdepOrder = 0, rho = 0, rand.gen = c("rnorm", "rpois"), noise.sd = .1)
{
  # Recycling necessary parameters
  ncat <- rep_len(ncat, p)
  fact <- !is.na(ncat)
  ffuns <- rep_len(ffuns, p)
  fbetas <- rep_len(fbetas, p)
  s <- rep_len(s, p)
  obetas <- rep_len(obetas, p + 1)
  extBetas <- rep_len(extBetas, p + 1)
  
  # Predictors X
  Smat <- matrix(rho, nrow = p, ncol = p)
  diag(Smat) <- 1
  Xsim <- MASS::mvrnorm(n, rep(0, p), Smat)
  Xsim <- as.data.frame(Xsim)
  if (any(fact)){
    Xsim[fact] <- Map(function(x, nc, s){
        quants <- quantile(x, c(seq(0, s, length.out = nc), 1))
        quants[1] <- quants[1] - 1
        quants[nc + 1] <- quants[nc + 1] + 1
        cut(x, quants, labels = 1:nc)
      }, 
      Xsim[fact], ncat[fact], s[fact]
    )
  }
  
  # To apply s
  stype <- match.arg(stype)
  if (stype == "quantile"){
    s <- mapply(quantile, Xsim, s)
  } 
  
  # Mean part
  Ypart <- Map(function(f, b, x) do.call(f, list(x = as.numeric(x), Betas = b)), 
    ffuns, fbetas, Xsim)
  obetas[-1] <- obetas[-1] / p
  linPred <- cbind(1, data.matrix(as.data.frame(Ypart))) %*% obetas
  
  # Detect extremes
  uni.extremes <- matrix(NA, nrow = n, ncol = p)
  uni.extremes[,!fact] <- mapply(">=", Xsim[!fact], s[!fact])
  if (any(fact)) uni.extremes[,fact] <- mapply("==", Xsim[fact], 
    sapply(Xsim[fact], function(x) max(as.numeric(x))))
  extremes <- apply(uni.extremes, 1, all)
  
  # Apply extreme function
  Xext <- matrix(NA, nrow = n, ncol = p)
  Xext[,!fact] <- mapply("-", Xsim[!fact], s[!fact])
  Xext[,!fact] <- mapply("/", as.data.frame(Xext)[,!fact, drop = F], 
    sapply(Xsim[!fact], max) - s[!fact])
  if (any(fact)){
    Xext[,fact] <- sapply(Xsim[fact], 
      function(x) as.numeric(x == max(as.numeric(x))))
  }
  extBetas[-1] <- extBetas[-1] / p
  extPred <- ifelse(extremes, cbind(1, Xext) %*% extBetas, 0)

  # Response simulation
  rand.gen <- match.arg(rand.gen)
  simpars <- list(model = list(ar = YdepOrder), n = n, rand.gen = get(rand.gen))
  if (rand.gen == "rnorm"){
    simpars$mean <- linPred + extPred
    simpars$sd <- noise.sd
  } else {
    simpars$lambda <- exp(linPred + extPred)
  }
  suppressWarnings(Ysim <- do.call(arima.sim, simpars))

  # Return simulations
  output <- list(Ysim = Ysim, Xsim = Xsim, linPred = linPred, 
    extPred = extPred, extremes = extremes)
  return(output)
}
