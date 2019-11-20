# Computes the F-score using observed and predicted indices of the event.
#
# The F-score is based on recall which is identical to sensivity,
#  i.e. the proportion of events predicted, and sensitivity which is the 
#  proportion of predicted events that are true events. 
#
# predicted   Vector of index of predicted events.
# observed    Vector of index of observed events.
# beta        Beta value of the F-score. If larger than 1, more importance
#  is given to recall while if lower, more importance is given to precision.
Fscore <- function(predicted, observed, beta = 1)
{
  nObs <- length(observed)
  nPred <- length(predicted)
  truePos <- sum(predicted %in% observed)
  recall <- truePos / nObs
  precision <- truePos / nPred
  Fscore <- (precision * recall) / ((beta^2 * precision) + recall)
  return(Fscore)
}

# Find local extrema and zero-crossings of a data series
#
# Return a list giving the indices of local minima, 
#   local maxima and the total number of extrema.
#
# x   The signal.
find.extrema <- function(x){
    n <- length(x)
    # when consecutive points have the same value, keep only one (the one on the middle)
    adoub <- which(diff(x) != 0)
    inddoub <- which(diff(adoub) != 1) + 1
    d <- adoub[inddoub] - adoub[inddoub-1]
    adoub[inddoub] <- adoub[inddoub] - floor(d/2)
    adoub <- c(adoub,n)
    x1 <- x[adoub]
    # Search the extrema among the remaining points
    d2x <- diff(sign(diff(x1)))
    indmin <- adoub[which(d2x > 0) + 1] 
    indmax <- adoub[which(d2x < 0) + 1]
    nextr <- c(length(indmin),length(indmax))
    names(nextr) <- c("nmin","nmax")
    return(list(indmin = indmin, indmax = indmax, nextrema = nextr))
}

# Extracts the thresholds from a party object (see partykit)
#
# tree    A party object as returned by functions lmtree or ctree
extractThresholds_party <- function(tree){
  dat <- tree$data
  yind <- attr(tree$terms, "response")
  xnames <- attr(tree$info$terms$partitioning, "term.labels")
  
  # Find the highest mean node
  nodeMean <- aggregate(dat[,yind], 
    by = list(node = predict(tree, type = "node")), 
    mean)
  nodeNum <- with(nodeMean, node[which.max(x)])
  
  # Extract thresholds
  nodepath <- partykit:::.list.rules.party(tree, nodeNum)
  inner_path <- unlist(strsplit(nodepath, " & "))
  rules <- strsplit(inner_path, "<|>=?")
  rules <- t(as.data.frame(rules))
  findthresh <- aggregate(as.numeric(rules[,2]), by = list(var = rules[,1]), 
    max)
  rownames(findthresh) <- findthresh[,1]
  thresholds <- findthresh[xnames, 2]
  
  return(list(thresholds = thresholds, bestnode = nodeNum))
}

# Extracts the thresholds from a gam object (see mgcv)
#
# object    A gam object as returned by the function gam
extractThresholds_gam <- function(object){
  # Compute lower 95% confidence bound
  exp_fun <- predict(object, type = "terms", se.fit = TRUE)
  ci_lo <- exp_fun$fit - 1.96 * exp_fun$se.fit
  
  # Sort matrices according to exposure
  xb <- object$model[,-attr(object$terms, "response")]
  ords <- apply(xb, 2, order)
  xb_ord <- Map("[", as.data.frame(xb), as.data.frame(ords))
  fit_ord <- Map("[", as.data.frame(exp_fun$fit), as.data.frame(ords))
  ci_ord <- Map("[", as.data.frame(ci_lo), as.data.frame(ords))
  
  # Determine minimum mortality temperature as the highest local minimum in the
  #   exposure response function
  mmt <- mapply(function(f, x){
    f <- round(f, digits = getOption("digits")) # Avoid small numerical
                                      # imprecisions to perturb MMT finding
    ind <- max(1, emdr:::find.extrema(f)$indmin)
    return(x[ind])
  }, fit_ord, xb_ord)
  above_mmt <- mapply(">", xb_ord, mmt) 
  
  # Estimate threshold as the lowest temperature with lo_ci > 0
  signif_fun <- sapply(ci_ord, ">", 0)
  thresh_ind <- apply(above_mmt & signif_fun, 2, function(x) min(which(x)))
  thresholds <- mapply("[", xb_ord, thresh_ind)
  
  return(thresholds)
}

# Detects the alarms from the data using the given thresholds
#
# thresholds    A vector of thresholds. Must be consistent with the number
#   of columns in x
# x   The indicators
# y   The response
extract_alarms <- function(thresholds, x, y){  
  x <- data.frame(x)
  uni.alb <- mapply(">=", x, thresholds)
  alarmlogi <- apply(uni.alb, 1, all, na.rm = T)
  alarms <- y[alarmlogi]
  names(alarms) <- which(alarmlogi)
  return(alarms)
}

# Transforms the line coordinate in user coordinates in a plot margins
#
# line    Line coordinates
# side    The margin side
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
    `1` = par('usr')[3] - line * y_off * lh,
    `2` = par('usr')[1] - line * x_off * lh,
    `3` = par('usr')[4] + line * y_off * lh,
    `4` = par('usr')[2] + line * x_off * lh,
    stop("side must be 1, 2, 3, or 4", call.=FALSE))
}

# Add a legend in the margins of a plot
#
# location   A character string giving the location of the legend. One of:
#     "topcenter", "topleft", "topright", "rightcenter", "righttop", 
#     "rightbottom", "leftcenter", "lefttop", "leftbottom", "bottomcenter", 
#     "bottomleft", "bottomright"
# ...    The usual parameters of the function legend
outerLegend <- function(location,...)
{
  if (length(dev.list())==0) stop("Aucun graphique ouvert") 
  olparams <- switch(tolower(location),
    topcenter = list(originalLocation = "top", direction = c(0,1)),
    topright = list(originalLocation = "topright", direction = c(0,1)),
    topleft = list(originalLocation = "topleft", direction = c(0,1)),
    rightcenter = list(originalLocation = "right", direction = c(1,0)),
    righttop = list(originalLocation = "topright", direction = c(1,0)),
    rightbottom = list(originalLocation = "bottomright", direction = c(1,0)),
    leftcenter = list(originalLocation = "left", direction = c(1,0)),
    lefttop = list(originalLocation = "topleft", direction = c(1,0)),
    leftbottom = list(originalLocation = "bottomleft", direction = c(1,0)),
    bottomcenter = list(originalLocation = "bottom", direction = c(0,1)),
    bottomright = list(originalLocation = "bottomright", direction = c(0,1)),
    bottomleft = list(originalLocation = "bottomleft", direction = c(0,1)),
    stop("Unknown 'location' parameter")
  )
  leg <- legend(olparams$originalLocation,plot=F,...) 
  insetx <- leg$rect$w/(par("usr")[2]-par("usr")[1])
  insety <- leg$rect$h/(par("usr")[4]-par("usr")[3])
  insetpar <- -c(insetx,insety) * olparams$direction
  legend(olparams$originalLocation,xpd=NA,inset=insetpar,...)
}