#########################################################
#
#            Parametric simulations
#         Part 0: Simulation functions
#                                                  
#########################################################

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
#' @param s Thresholds to be found. 
#' @param stype If "quantile" s corresponds to the quanile of X. If not, directly
#'   directly corresponds to the threshold.
#' @param extBetas Parameters for the linear function of extremes.
#' @param XdepOrder AR order for temporal dependence of indicators.
#' @param YdepOrder AR ordre for temporal dependence of response.
#' @param noise.sd Standard deviation of the random part of response, relative
#'    to the standard deviation of the deterministic part.
generate.data <- function(B = 5000, n = 5000, p = 2, Lsim = 1, 
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
  
  Xext <- mapply("-", as.data.frame(Xsim), s)
  Xext <- mapply("/", as.data.frame(Xext), apply(Xsim, 2, max) - s)
  extPred <- ifelse(extremes, cbind(1, Xext) %*% extBetas, 0)

  # Response simulation
  suppressWarnings(Ysim <- replicate(B, scale(linPred + extPred) + 
    arima.sim(list(ar = YdepOrder), n, rnorm, sd = noise.sd), simplify = TRUE))
  
  output <- list(Ysim = Ysim, Xsim = Xsim, linPred = linPred, extPred = extPred,
    fTrue = fTrue, overall = overall, parameters = params, s = s)
  return(output)
}
