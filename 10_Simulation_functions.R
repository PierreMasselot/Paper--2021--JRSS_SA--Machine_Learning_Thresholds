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
#' @param stype If "quantile" s corresponds to the quantile of X. If not, 
#'   directly corresponds to the threshold.
#' @param extBetas Parameters for the linear function of extremes.
#' @param XdepOrder AR order for temporal dependence of indicators.
#' @param YdepOrder AR ordre for temporal dependence of response.
#' @param noise.sd Standard deviation of the random part of response, relative
#'    to the standard deviation of the deterministic part.
generate.data <- function(n = 5000, p = 2, ncat = NA, 
  obetas = 1, ffuns = "fconstant", fbetas = list(), 
  s = .8, stype = c("quantile", "absolute"), extBetas = 1,
  YdepOrder = 0, rho = 0, noise.sd = .2)
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
    Xsim[!fact] <- lapply(Xsim[!fact], function(x) rank(x) / n)
  } else {
    Xsim[!fact] <- lapply(Xsim[!fact], scale)
  }
  
  # Mean part
  Ypart <- Map(function(f, b, x) do.call(f, list(x = as.numeric(x), Betas = b)), 
    ffuns, fbetas, Xsim)
  linPred <- cbind(1, scale(data.matrix(as.data.frame(Ypart)))) %*% obetas
  linPred <- scale(linPred)
  
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
  extPred <- ifelse(extremes, cbind(1, Xext) %*% extBetas, 0)

  # Response simulation
  suppressWarnings(Ysim <- scale(linPred + extPred) + 
    arima.sim(list(ar = YdepOrder), n, rnorm, sd = noise.sd))
  
  output <- list(Ysim = Ysim, Xsim = Xsim, linPred = linPred, 
    extPred = extPred, extremes = extremes)
  return(output)
}
