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

MOB.apply <- function(yb, xb, zb = NULL, ...){
  # Data and formula
  datab <- data.frame(y = yb, x = xb)
  if (!is.null(zb)) datab <- data.frame(datab, z = zb)
  datab <- na.omit(datab)
  xinds <- grep("x\\.?", colnames(datab))
  zinds <- grep("z\\.?", colnames(datab))
  
  # Grow the tree (glmtree is a wrapper for mob with GLM)
  modform <- paste(colnames(datab)[zinds], collapse = " + ")
  splitform <- paste(colnames(datab)[xinds], collapse = " + ")
  form <- sprintf("y ~ %s %s %s", modform, ifelse(is.null(zb), "", "|"), 
    splitform)   
  treeb <- glmtree(as.formula(form), data = datab, ...)
  
  # Extract thresholds
  extractThresholds_party(treeb)$thresholds
}


MARS.apply <- function(yb, xb, zb = NULL, ...){
  datab <- data.frame(y = yb, x = xb)
  if (!is.null(zb)) datab <- data.frame(datab, z = zb)
  datab <- na.omit(datab)
  xinds <- grep("x\\.?", colnames(datab))
  
  # Apply MARS
  marsb <- earth(y ~ ., data = datab, ...)
  
  # Thresholds
  apply(marsb$cuts[,xinds - 1, drop = FALSE], 2, max)
}


PRIM.apply <- function(yb, xb, zb = NULL, RRind = 1:ncol(zb), 
  family = "gaussian", ...){
  keep <- complete.cases(yb)
  xb <- xb[keep,, drop = FALSE]
  yb <- na.omit(yb)
  n <- length(yb)
  if (is.null(zb)){
    obj.fun <- mean
  } else {
    zb <- zb[keep,, drop = FALSE]
    obj.fun <- function(y, x, inbox){
      y <- y[inbox]
      dat <- data.frame(y, zb[inbox,])
      fit <- glm(y ~ ., data = dat, family = family)
      pred <- coef(fit)[RRind + 1]
      return(mean(pred))
    }
  }
  
  # Peeling the box
  peelb <- peeling(yb, xb, obj.fun = obj.fun, ...)
  
  # Extract box : box following the biggest increase in mean
  chosenb <- jump.prim(peelb)  
  boxb <- pasting(peelb, support = chosenb$final.box$support)
  
  sapply(extract.box(boxb)$limits[[1]], "[", 1)    
}


AIM.apply <- function(yb, xb, zb = NULL, numcut = 3, ...){  
  keep <- complete.cases(yb)
  yb <- na.omit(yb)
  if (!is.null(zb)){
    xb <- data.frame(x = xb, zb)
  } else {
    xb <- data.frame(x = xb)
  }
  xb <- xb[keep,, drop = FALSE]
  n <- length(yb)
  xinds <- grep("x\\.?", colnames(xb))
  p <- length(xinds)
  
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
    return(rep(NA_real_, p))
  }
  
  aimb <- lm.main(xb, yb, nsteps = cvb$kmax,
    maxnumcut = numcut, ...)
  
  # Extract thresholds
  rules <- aimb$res[[cvb$kmax]]
  maxrules <- aggregate(rules[,2], by = list(var = rules[,1]), max)
  thresholds <- rep(NA_real_, p)
  thresholds[xinds %in% maxrules$var] <- maxrules$x[maxrules$var %in% xinds]
  
  return(thresholds)
}


GAM.apply <- function(yb, xb, zb = NULL, ...){
  datab <- data.frame(y = yb, x = xb)
  if (!is.null(zb)) datab <- cbind(datab, zb)
  datab <- na.omit(datab)
  xinds <- grep("x\\.?", colnames(datab))
  
  # Fitting the GAM
  form_rhs <- paste(sprintf("s(%s)", colnames(datab)[-1]), collapse = " + ")
  gam_form <- sprintf("y ~ %s", form_rhs)
  gam_res <- gam(as.formula(gam_form), data = datab, ...)
  
  extractThresholds_gam(gam_res)[xinds - 1]
}

DLNM.apply <- function(yb, xb, zb = NULL, crosspars = vector("list", ncol(xb)), gampars = list()){
  # Prepare crossbases
  cblist <- vector("list", ncol(xb))
  for (i in seq_len(ncol(xb))){
    crosspars[[i]]$x <- xb[,i]
    cblist[[i]] <- do.call(crossbasis, crosspars[[i]])
  }
  names(cblist) <- sprintf("cb%i", 1:length(cblist))
  
  # Prepare penalties
  bfuns <- sapply(cblist, function(x) attr(x, "argvar")$fun)
  if (any(bfuns %in% c("ps", "cr"))){
    penlist <- lapply(cblist, cbPen)
    gampars$paraPen <- penlist
  }
  
  # Fit model
  datab <- c(list(y = yb), cblist)
  form <- sprintf("y ~ %s", paste(names(cblist), collapse = " + "))
  if (!is.null(zb)) form <- sprintf("%s + zb", form)
  gampars$formula <- as.formula(form)
  gampars$data <- datab
  res <- do.call(gam, gampars)
  
  # Obtain overall relationship
  gamcoefs <- res$coefficients
  thresholds <- vector("numeric", ncol(xb))
  for (i in seq_len(ncol(xb))){
    inds <- grep(names(cblist)[i], names(coef(res)))
    # First estimate of the association
    cp <- crosspred(cblist[[i]], coef = gamcoefs[inds], vcov = vcov(res)[inds, inds, drop = F], 
      by = .001)
    thresholds[i] <- min(cp$predvar[cp$alllow > 0])
    # Minimum mortality temperature
    mmtind <- max(1, emdr:::find.extrema(cp$allfit)$indmin)
    mmt <- cp$predvar[mmtind]
    # Find the lowest value (above mmt) such that alllow is above 0
    abovemmt <- cp$predvar > mmt
    thresholds[i] <- min(cp$predvar[abovemmt & cp$alllow > 0])
  }
  
  return(thresholds)
}


#chngpt.apply <- function(yb, xb, zb = NULL, ...){
#  datab <- data.frame(y = yb, x = xb)
#  if (!is.null(zb)) datab <- cbind(datab, zb)
#  datab <- na.omit(datab)
#  xinds <- grep("x\\.", colnames(datab))
#  
#  # Prepare formulas
#  f1 <- sprintf("y ~ %s", ifelse(is.null(zb), "1", paste(colnames(zb), collapse = "+")))
#  f2 <- sprintf("~ %s", paste(colnames(datab)[xinds], collapse = " + "))
#  
#  # Fit threshold regression model
#  thresh_res <- chngptm(as.formula(f1), as.formula(f2), data = datab, 
#    family = "gaussian", type = "stegmented")
#}

seg.apply <- function(yb, xb, zb = NULL, glmpars = list(), segpars = list())
{
  # Data
  datab <- data.frame(y = yb, x = xb)
  if (!is.null(zb)) datab <- data.frame(datab, z = zb)
  datab <- na.omit(datab)
  xinds <- grep("x\\.?", colnames(datab))
  zinds <- grep("z\\.?", colnames(datab))
  
  # Basic model
  modrhs <- ifelse(is.null(zb), "1", paste(names(datab)[zinds], collapse = " + "))
  modform <- sprintf("y ~ %s", modrhs)
  glmpars$formula <- as.formula(modform)
  glmpars$data <- datab
  mod <- do.call(glm, glmpars)
  
  # Segmented
  segpars$obj <- mod
  segform <- sprintf("~ %s", paste(names(datab)[xinds], collapse = " + "))
  segpars$seg.Z <- as.formula(segform)
  resb <- do.call(segmented, segpars)
  
  # Extract thresholds
  return(resb$psi[,2])
}