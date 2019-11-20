####################################################
#
#             Threshold finding methods
#
####################################################

# Wrapper functions for the threshold estimation methods. All functions
#   contain the same parameters:
# yb    The response vector
# xb    The indicator matrix
# ...   Parameters to be passed to the main method (e.g. lmtree, earth, ...)

MOB.apply <- function(yb, xb, ...){
  # Data and formula
  datab <- data.frame(y = yb, x = xb)
  xinds <- grep("x\\.", colnames(datab))
  
  # Grow the tree (lmtree is a wrapper for mob with linear regression)
  form <- sprintf("y ~ %s", paste(colnames(datab)[-1], collapse = " + "))   
  treeb <- lmtree(as.formula(form), data = datab, ...)
  
  # Extract thresholds
  thresholds <- extractThresholds_party(treeb)
  
  # Alarms
  alarms <- extract_alarms(thresholds, xb, yb)
  
  return(list(thresholds = thresholds$thresholds, alarms = alarms))
}


MARS.apply <- function(yb, xb, ...){
  datab <- data.frame(y = yb, x = xb)
  xinds <- grep("x\\.", colnames(datab)) 
  
  # Apply MARS
  marsb <- earth(y ~ ., data = na.omit(datab), ...)
  
  # Thresholds
  thresholds <- apply(marsb$cuts[,xinds - 1], 2, max)
  
  # Extract alarms
  alarms <- extract_alarms(thresholds, xb, yb)
  
  list(thresholds = thresholds, alarms = alarms)
}


PRIM.apply <- function(yb, xb, ...){
  keep <- complete.cases(yb)
  xb <- xb[keep,]
  yb <- na.omit(yb)
  n <- length(yb)
  
  # Peeling the box
  peelb <- peeling(yb, xb, ...)
  
  # Extract box : box following the biggest increase in mean
  chosenb <- jump.prim(peelb)
  
  boxb <- pasting(peelb, support = chosenb$final.box$support)
  
  # Extract alarms
  alarmb <- predict(boxb)
  alarms <- yb[alarmb$inbox]
  names(alarms) <- which(alarmb$inbox)

  thresholds <- sapply(extract.box(boxb)$limits[[1]], "[", 1)    
  
  list(thresholds = thresholds, alarms = alarms)
}


AIM.apply <- function(yb, xb, numcut = 3, ...){
  xb <- data.matrix(xb)
  xb <- xb[!is.na(yb),]
  yb <- na.omit(yb)
  p <- ncol(xb)
  n <- length(yb)
  
  # Fit the model
  while (numcut > 0){
    cvb <- try(cv.lm.main(xb, yb, nsteps = numcut * p, 
      maxnumcut = numcut, ...), silent = TRUE)
    if (inherits(cvb, "try-error")){
      numcut <- numcut - 1
    } else {
      break
    }
  }
  if (inherits(cvb, "try-error")){
    return(list(thresholds = rep(NA_real_, p), alarms = NULL))
  }
  
  aimb <- lm.main(xb, yb, nsteps = cvb$kmax,
    maxnumcut = numcut, ...)
  
  # Extract thresholds
  rules <- aimb$res[[cvb$kmax]]
  maxrules <- aggregate(rules[,2], by = list(var = rules[,1]), max)
  thresholds <- rep(NA_real_, p)
  thresholds[maxrules$var] <- maxrules$x
  
  # Alarms
  alarms <- extract_alarms(thresholds, xb, yb)
  
  return(list(thresholds = thresholds, alarms = alarms))
}


GAM.apply <- function(yb, xb){
  datab <- data.frame(Y = yb, x = xb)
  
  # Fitting the GAM
  form_rhs <- paste(sprintf("s(%s)", colnames(datab)[-1]), collapse = " + ")
  gam_form <- sprintf("Y ~ %s", form_rhs)
  gam_res <- gam(as.formula(gam_form), data = datab)
  
  thresholds <- extractThresholds_gam(gam_res)
  
  # Extract alarms
  alarms <- extract_alarms(thresholds, xb, yb)
  
  list(thresholds = thresholds, alarms = alarms)
}


