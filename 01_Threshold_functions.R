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
  xinds <- 1:ncol(xb) + 1
  if (!is.null(zb)) {
    datab <- data.frame(datab, z = zb)
    zinds <- 1:ncol(zb) + ncol(xb) + 1
  }
  datab <- na.omit(datab)
  
  # Grow the tree (glmtree is a wrapper for mob with GLM)
  modform <- paste(colnames(datab)[zinds], collapse = " + ")
  splitform <- paste(colnames(datab)[xinds], collapse = " + ")
  form <- sprintf("y ~ %s %s %s", modform, ifelse(is.null(zb), "", "|"), 
    splitform)   
  treeb <- glmtree(as.formula(form), data = datab, ...)
  
  # Extract thresholds
  thresholds <- extractThresholds_party(treeb)
  
  # Extract alerts
  alerts <- which(predict(treeb, type = "node") == thresholds$bestnode)
  
  list(thresholds = thresholds$thresholds, alerts = alerts)
}


MARS.apply <- function(yb, xb, zb = NULL, endspan, ...){
  datab <- data.frame(y = yb, x = xb)
  xinds <- 1:ncol(xb)
  if (!is.null(zb)) {
    datab <- data.frame(datab, z = zb)
  }
  datab <- na.omit(datab)
  
  # Apply MARS
  marsb <- earth(y ~ ., data = datab, endspan = endspan, ...)
  
  # Thresholds
  cutsbyvar <- lapply(as.data.frame(marsb$cuts[,xinds, drop = FALSE]), 
    function(x) sort(unique(x[x != 0])))
  nocut <- sapply(cutsbyvar, length) == 0
  if (any(nocut)){
    cutsbyvar[nocut] <- lapply(datab[xinds[nocut] + 1], min)
  }
  
  # try different threshold combination
  allthreshs <- expand.grid(cutsbyvar)
  nalerts <- apply(allthreshs, 1, 
    function(x) sum(apply(mapply(">=", as.data.frame(xb), x), 1, all))
  )
  meanalerts <- apply(allthreshs, 1, 
    function(x) mean(yb[apply(mapply(">=", as.data.frame(xb), x), 1, all)])
  )
  
  # Take the threshold with highest response among those with at least endspan
  #   observations
  ind <- which(rank(meanalerts) == max(rank(meanalerts)[nalerts >= endspan]))
  ind <- ind[length(ind)]
  if (length(ind) == 0) ind <- which.max(nalerts)
  thresholds <- allthreshs[ind,]
  
  # Alerts
  above <- mapply(">=", as.data.frame(xb), thresholds)
  alerts <- which(apply(above, 1, all))
  
  thresholds[nocut] <- NA
  
  list(thresholds = unlist(thresholds), alerts = alerts)
}


PRIM.apply <- function(yb, xb, zb = NULL, RRind = 1:ncol(zb), 
  family = "gaussian", ...)
{
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
  
  # Extract thresholds
  thresholds <- sapply(extract.box(boxb)$limits[[1]], "[", 1) 
  
  # Extract predicted extremes
  extremes <- predict(boxb)
  
  list(thresholds = thresholds, alerts = which(extremes$inbox))
}


AIM.apply <- function(yb, xb, zb = NULL, numcut = 3, mincut, ...){  
  keep <- complete.cases(yb)
  yb <- na.omit(yb)
  xinds <- 1:ncol(as.matrix(xb))
  if (!is.null(zb)){
    xb <- data.frame(x = xb, zb)
  } else {
    xb <- data.frame(x = xb)
  }
  xb <- xb[keep,, drop = FALSE]
  n <- length(yb)
  p <- length(xinds)
  # To avoid error in lm.main
  if (p == 1) xb <- cbind(xb, 1)
  
  # Fit the model
  while (numcut > 0){
    cvb <- try(cv.lm.main(xb, yb, nsteps = numcut * p, 
      maxnumcut = numcut, mincut = mincut, ...), silent = TRUE)
    if (inherits(cvb, "try-error")){
      numcut <- numcut - 1
    } else {
      break
    }
  }
  if (inherits(cvb, "try-error")){
    return(list(thresholds = rep(NA_real_, p),
      alerts = seq_along(yb)))
  }
  
  aimb <- lm.main(xb, yb, nsteps = cvb$kmax,
    maxnumcut = numcut, mincut = mincut, ...)
  
  # Extract thresholds
  rules <- aimb$res[[cvb$kmax]]
  rules <- rbind(rules, cbind(xinds, apply(xb[,xinds, drop = F], 2, min), 1))
  
  # Thresholds
  rulesbyvar <- split(rules[,2], rules[,1])[xinds]
  rulesbyvar <- lapply(rulesbyvar, sort)
  
  # try different threshold combination
  allthreshs <- expand.grid(rulesbyvar)
  nalerts <- apply(allthreshs, 1, 
    function(x) sum(apply(mapply(">=", as.data.frame(xb)[,xinds, drop = F], x), 
      1, all))
  )
  meanalerts <- apply(allthreshs, 1, 
    function(x) mean(yb[apply(mapply(">=", as.data.frame(xb)[,xinds, drop = F], 
      x), 1, all)])
  )
  
  # Take the threshold with highest response among those with at least endspan
  #   observations
  ind <- which(rank(meanalerts) == max(rank(meanalerts)[nalerts >= mincut]))
  ind <- ind[length(ind)]
  if (length(ind) == 0) ind <- which.max(nalerts)
  thresholds <- allthreshs[ind,]
  
  # Alerts
  above <- mapply(">=", as.data.frame(xb[,xinds, drop = F]), thresholds)
  alerts <- which(apply(above, 1, all))
  
  # Remove threshold
  thresholds[thresholds == apply(xb[,xinds, drop = F], 2, min)] <- NA
  
  list(thresholds = unlist(thresholds), alerts = alerts)
}


GAM.apply <- function(yb, xb, zb = NULL, ...){
  datab <- data.frame(y = yb, x = xb)
  xinds <- 2:ncol(datab)
  if (!is.null(zb)) datab <- cbind(datab, zb)
  datab <- na.omit(datab)
  
  # Fitting the GAM
  form_rhs <- paste(sprintf("s(%s)", colnames(datab)[-1]), collapse = " + ")
  gam_form <- sprintf("y ~ %s", form_rhs)
  gam_res <- gam(as.formula(gam_form), data = datab, ...)
  
  # extract thresholds
  thresholds <- extractThresholds_gam(gam_res)[xinds - 1]
  
  # Alerts
  above <- mapply(">=", as.data.frame(xb), thresholds)
  alerts <- which(apply(above, 1, all))
  
  list(thresholds = thresholds, alerts = alerts)
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
    inds <- grep(paste0(names(cblist)[i], "v"), names(coef(res)))
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
  xinds <- 2:ncol(datab)
  names(datab)[-1] <- sprintf("X%i", xinds - 1)
  if (!is.null(zb)) {
    datab <- data.frame(datab, z = zb)
    zinds <- 1:ncol(zb) + max(xinds)
    names(datab)[zinds] <- sprintf("Z%i", 1:ncol(zb))
  }
  datab <- na.omit(datab)
  
  
  # Basic model
  modrhs <- ifelse(is.null(zb), "1", paste(names(datab)[zinds], collapse = " + "))
  modform <- sprintf("y ~ %s", modrhs)
  glmpars$formula <- as.formula(modform)
  glmpars$data <- datab
  mod <- do.call(glm, glmpars)
  
  # Segmented
  if (!is.null(segpars$npsi)){
    segpars$npsi <- rep(segpars$npsi, length(xinds))
    names(segpars$npsi) <- names(datab)[xinds]
  } 
  if (!is.null(segpars$psi)){
    segpars$psi <- rep_len(segpars$psi, length(xinds))
    names(segpars$psi) <- names(datab)[xinds]
  } 
  segpars$obj <- mod
  segform <- sprintf("~ %s", paste(names(datab)[xinds], collapse = " + "))
  segpars$seg.Z <- as.formula(segform)
  resb <- try(do.call(segmented, segpars))
  
  # Extract thresholds
  if (inherits(resb, "try-error") || is.null(resb$psi)){    
    return(list(thresholds = rep(NA_real_, length(xinds)), 
      alerts = seq_along(yb)))
  } else {
    # Extract thresholds
    psivars <- sapply(strsplit(rownames(resb$psi), "\\."), "[", 2) 
    psires <- tapply(resb$psi[,2], psivars, max)
    names(psires) <- NULL
    
    # Alerts
    above <- mapply(">=", as.data.frame(xb), psires)
    alerts <- which(apply(above, 1, all))
    
    return(list(thresholds = psires, alerts = alerts))
  }  
}