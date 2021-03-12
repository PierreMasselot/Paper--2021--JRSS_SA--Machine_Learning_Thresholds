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
  if (nodeNum == 1) return(list(
    thresholds = rep(NA_real_, length(xnames)), 
    bestnode = nodeNum))
  
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
  return(list(indmin=indmin,indmax=indmax,nextrema=nextr))
}

# Extracts the thresholds from a gam object (see mgcv)
#
# object    A gam object as returned by the function gam
extractThresholds_gam <- function(object){
  # Compute lower 95% confidence bound
  exp_fun <- predict(object, type = "terms", se.fit = TRUE)
  ci_lo <- exp_fun$fit - 1.96 * exp_fun$se.fit
  
  # Sort matrices according to exposure
  xb <- object$model[,-attr(object$terms, "response"), drop = FALSE]
  ords <- apply(xb, 2, order)
  xb_ord <- Map("[", as.data.frame(xb), as.data.frame(ords))
  fit_ord <- Map("[", as.data.frame(exp_fun$fit), as.data.frame(ords))
  ci_ord <- Map("[", as.data.frame(ci_lo), as.data.frame(ords))
  
  # Determine minimum mortality temperature as the highest local minimum in the
  #   exposure response function
  mmt <- mapply(function(f, x){
    f <- round(f, digits = getOption("digits")) # Avoid small numerical
                                      # imprecisions to perturb MMT finding
    ind <- max(1, find.extrema(f)$indmin)
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

axis.intervals <- function(side=1, ticks = axTicks(side), atLabels = NULL, labels = 1:length(atLabels), ...)
# side : integer between 1 and 4, the side at which to draw the interval.  The axis is placed as follows: 1=below, 2=left, 3=above and 4=right ;
# ticks : position of the ticks delimiting the intervals. If a single number is provided, this indicates the number of intervals; if a vector is provided, this indicates the position of each ticks ; 
# atLabels : the position of the labels for each interval. If not of length nIntervals - 1, the vector is either recycled, either cut. The default is to the center of the interval ;
# labels : the labels of each interval ; 
# ... : other graphical parameters to pass to the 'axis' function
{
    stopifnot((side <- as.integer(side)) %in% 1:4)
    if (length(ticks) == 1){
       is.x <- side%%2 == 1
       usr <- par("usr")[1:2 + 2*!is.x]
       XY <- function(ch) paste0(if (is.x) "x" else "y", ch)
       axs <- par(XY("axs"))
       if (axs == "r"){
          per4 <- 4*diff(usr)/108
          usr[1] <- usr[1] + per4
          usr[2] <- usr[2] - per4
       }
       nIntervals <- floor(ticks)
       ticks <- seq(usr[1],usr[2],length.out = nIntervals + 1)
    } else {
       nIntervals <- length(ticks) - 1
    }    
    axis(side, at = ticks, labels = FALSE, ...)
    if (is.null(atLabels)){
       atLabels <- (ticks[-1] + ticks[1:nIntervals]) / 2
    } else {
       if (length(atLabels) != nIntervals) atLabels <- rep_len(atLabels,nIntervals)
    }
    if (length(labels) != nIntervals) labels <- rep_len(labels,nIntervals)
    axis(side, at = atLabels, labels = labels, tick = FALSE, ...)
}

node_bivplotmod <- function (mobobj, which = NULL, id = TRUE, pop = TRUE, pointcol = "black", 
  pointcex = 0.5, boxcol = "black", boxwidth = 0.5, boxfill = "lightgray", 
  bg = "white", fitmean = TRUE, linecol = "red", lwd = 1, pch = 1,
  cdplot = FALSE, fivenum = TRUE, breaks = NULL, ylines = NULL, 
  xlab = FALSE, ylab = FALSE, margins = rep(1.5, 4), mainlab = NULL, 
  ...) 
{
  mf <- model.frame(mobobj)
  y <- Formula::model.part(mobobj$info$Formula, mf, lhs = 1L, 
    rhs = 0L)
  if (isTRUE(ylab)) 
    ylab <- names(y)
  if (identical(ylab, FALSE)) 
    ylab <- ""
  if (is.null(ylines)) 
    ylines <- ifelse(identical(ylab, ""), 0, 2)
  y <- y[[1L]]
  X <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, 
    rhs = 1L)
  fitted <- mobobj$fitted[["(fitted)"]]
  if (inherits(X, "try-error")) {
    rval <- switch(class(y)[1L], Surv = node_surv(mobobj, 
      id = id, mainlab = mainlab, ...), factor = node_barplot(mobobj, 
        id = id, mainlab = mainlab, ...), ordered = node_barplot(mobobj, 
          id = id, mainlab = mainlab, ...), node_boxplot(mobobj, 
            ...))
    return(rval)
  }
  if (is.factor(y)) 
    y <- factor(y, levels = rev(levels(y)))
  if (is.null(which)) 
    which <- 1L:NCOL(X)
  X <- X[, which, drop = FALSE]
  k <- NCOL(X)
  xlab <- if (!identical(xlab, FALSE)) {
    if (isTRUE(xlab)) 
      colnames(X)
    else rep(xlab, length.out = k)
  }
  else rep("", k)
  if (is.factor(y)) {
    if (!requireNamespace("vcd")) 
      stop(sprintf("Package %s is required for spine/CD plots", 
        sQuote("vcd")))
    if (cdplot) {
      num_fun <- function(x, y, yfit, i, name, ...) {
        vcd::cd_plot(x, y, xlab = xlab[i], ylab = ylab, 
          name = name, newpage = FALSE, margins = margins, 
          pop = FALSE, ...)
        if (fitmean) {
          grid.lines(x, yfit, default.units = "native", 
            gp = gpar(col = linecol, lwd = lwd))
          if (pop) 
            popViewport()
          else upViewport()
        }
        else {
          if (pop) 
            popViewport()
          else upViewport()
        }
      }
    }
    else {
      xscale <- if (is.null(breaks)) {
        if (fivenum) 
          lapply(X, function(z) {
            if (is.factor(z)) 
              1
            else fivenum(z)
          })
        else lapply(X, function(z) {
          if (is.factor(z)) 
            1
          else hist(z, plot = FALSE)$breaks
        })
      }
      else {
        if (is.list(breaks)) 
          breaks
        else list(breaks)
      }
      num_fun <- function(x, y, yfit, i, name, ...) {
        vcd::spine(x, y, xlab = xlab[i], ylab = ylab, 
          name = name, newpage = FALSE, margins = margins, 
          pop = FALSE, breaks = xscale[[i]], ...)
        if (fitmean) {
          xaux <- cut(x, breaks = xscale[[i]], include.lowest = TRUE)
          yfit <- unlist(tapply(yfit, xaux, mean))
          xaux <- prop.table(table(xaux))
          xaux <- cumsum(xaux) - xaux/2
          grid.lines(xaux, yfit, default.units = "native", 
            gp = gpar(col = linecol, lwd = lwd))
          grid.points(xaux, yfit, default.units = "native", 
            gp = gpar(col = linecol, cex = pointcex), 
            pch = pch)
          if (pop) 
            popViewport()
          else upViewport()
        }
        else {
          if (pop) 
            popViewport()
          else upViewport()
        }
      }
    }
    cat_fun <- function(x, y, yfit, i, name, ...) {
      vcd::spine(x, y, xlab = xlab[i], ylab = ylab, name = name, 
        newpage = FALSE, margins = margins, pop = FALSE, 
        ...)
      if (fitmean) {
        yfit <- unlist(tapply(yfit, x, mean))
        xaux <- prop.table(table(x))
        xaux <- cumsum(xaux + 0.02) - xaux/2 - 0.02
        grid.lines(xaux, yfit, default.units = "native", 
          gp = gpar(col = linecol, lwd = lwd))
        grid.points(xaux, yfit, default.units = "native", 
          gp = gpar(col = linecol, cex = pointcex), pch = pch)
        if (pop) 
          popViewport()
        else upViewport()
      }
      else {
        if (pop) 
          popViewport()
        else upViewport()
      }
    }
  }
  else {
    xscale <- sapply(X, function(z) {
      if (is.factor(z)) 
        c(1, length(levels(z)))
      else range(z)
    })
    yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
    num_fun <- function(x, y, yfit, i, name, ...) {
      xscale[, i] <- xscale[, i] + c(-0.1, 0.1) * diff(xscale[, 
        i])
      pushViewport(plotViewport(margins = margins, name = name, 
        yscale = yscale, xscale = xscale[, i]))
      grid.points(x, y, gp = gpar(col = pointcol, cex = pointcex), pch = pch)
      if (fitmean) {
        grid.lines(x, yfit, default.units = "native", 
          gp = gpar(col = linecol, lwd = lwd))
      }
      grid.xaxis(at = c(ceiling(xscale[1L, i] * 10), floor(xscale[2L, 
        i] * 10))/10)
      grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))
      grid.rect(gp = gpar(fill = "transparent"))
      if (ylab != "") 
        grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2.5, 
          "lines"), rot = 90)
      if (xlab[i] != "") 
        grid.text(xlab[i], x = unit(0.5, "npc"), 
          y = unit(-2, "lines"))
      if (pop) 
        popViewport()
      else upViewport()
    }
    cat_fun <- function(x, y, yfit, i, name, ...) {
      xlev <- levels(x)
      pushViewport(plotViewport(margins = margins, name = name, 
        yscale = yscale, xscale = c(0.3, xscale[2L, i] + 
            0.7)))
      for (j in seq_along(xlev)) {
        by <- boxplot(y[x == xlev[j]], plot = FALSE)
        xl <- j - boxwidth/4
        xr <- j + boxwidth/4
        grid.lines(unit(c(xl, xr), "native"), unit(by$stats[1L], 
          "native"), gp = gpar(col = boxcol, lwd = lwd))
        grid.lines(unit(j, "native"), unit(by$stats[1L:2L], 
          "native"), gp = gpar(col = boxcol, lty = 2, lwd = lwd))
        grid.rect(unit(j, "native"), unit(by$stats[2L], 
          "native"), width = unit(boxwidth, "native"), 
          height = unit(diff(by$stats[2:3]), "native"), 
          just = c("center", "bottom"), gp = gpar(col = boxcol, 
            fill = boxfill))
        grid.rect(unit(j, "native"), unit(by$stats[3L], 
          "native"), width = unit(boxwidth, "native"), 
          height = unit(diff(by$stats[3L:4L]), "native"), 
          just = c("center", "bottom"), gp = gpar(col = boxcol, 
            fill = boxfill))
        grid.lines(unit(j, "native"), unit(by$stats[4L:5L], 
          "native"), gp = gpar(col = boxcol, lty = 2, lwd = lwd))
        grid.lines(unit(c(xl, xr), "native"), unit(by$stats[5L], 
          "native"), gp = gpar(col = boxcol, lwd = lwd))
        n <- length(by$out)
        if (n > 0L) {
          grid.points(unit(rep.int(j, n), "native"), 
            unit(by$out, "native"), size = unit(0.5, 
              "char"), gp = gpar(col = boxcol), pch = pch)
        }
      }
      if (fitmean) {
        yfit <- unlist(tapply(yfit, x, mean))
        grid.lines(seq_along(xlev), yfit, default.units = "native", 
          gp = gpar(col = linecol, lwd = lwd))
        grid.points(seq_along(xlev), yfit, default.units = "native", 
          gp = gpar(col = linecol, cex = pointcex), pch = pch)
      }
      grid.rect(gp = gpar(fill = "transparent"))
      grid.xaxis(at = 1L:length(xlev), label = xlev)
      grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))
      if (ylab != "") 
        grid.text(ylab, y = unit(0.5, "npc"), x = unit(-3, 
          "lines"), rot = 90)
      if (xlab[i] != "") 
        grid.text(xlab[i], x = unit(0.5, "npc"), 
          y = unit(-2, "lines"))
      if (pop) 
        popViewport()
      else upViewport()
    }
  }
  rval <- function(node) {
    nid <- id_node(node)
    ix <- fitted %in% nodeids(mobobj, from = nid, terminal = TRUE)
    y <- y[ix]
    top_vp <- viewport(layout = grid.layout(nrow = k, ncol = 2, 
      widths = unit(c(ylines, 1), c("lines", "null")), 
      heights = unit(k, "null")), width = unit(1, 
        "npc"), height = unit(1, "npc") - unit(2, 
          "lines"), name = paste("node_mob", nid, 
            sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if (id) {
        function(id, nobs) sprintf("Node %s (n = %s)", 
          id, nobs)
      }
      else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, info_node(node)$nobs)
    }
    grid.text(mainlab, y = unit(1, "npc") - unit(0.75, 
      "lines"))
    popViewport()
    for (i in 1L:k) {
      xi <- X[ix, i]
      o <- order(xi)
      yi <- y[o]
      xi <- xi[o]
      yfit <- if (is.null(node$info$object)) {
        fitted(refit.modelparty(mobobj, node = nid))[o]
      }
      else {
        modmat <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, 
          rhs = 1L)[ix,]
        means <- apply(modmat, 2, mean)
        modmatnode <- matrix(means, nrow = sum(ix), ncol = length(means), 
          byrow = T)
        modmatnode[,grep(colnames(modmat)[i], names(means))] <- modmat[[i]]
        exp(cbind(1, modmatnode) %*% node$info$coefficients)[o]
      }
      plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = i)
      pushViewport(plot_vpi)
      if (is.factor(xi)) 
        cat_fun(xi, yi, yfit, i, paste("node_mob", 
          nid, "-", i, sep = ""), ...)
      else num_fun(xi, yi, yfit, i, paste("node_mob", 
        nid, "-", i, sep = ""), ...)
      if (pop) 
        popViewport()
      else upViewport()
    }
    if (pop) 
      popViewport()
    else upViewport()
  }
  return(rval)
}
class(node_bivplotmod) <- "grapcon_generator"
