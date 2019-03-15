####################################################
#
#             Threshold finding methods
#
####################################################

#' CART 
CART.apply <- function(yb, xb){
  datab <- data.frame(Y = yb, xb)
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
  uni.alb[is.na(uni.alb)] <- TRUE
  alarmb <- apply(uni.alb, 1, all)
  alarms <- datab[alarmb,1]
  names(alarms) <- which(alarmb)
   
  list(thresholds = boxb$box[1,], alarms = alarms)
}

#' MARS
#'
#' @param p Interaction degree
MARS.apply <- function(yb, xb, p = 2){
  datab <- data.frame(Y = yb, xb)
  # Apply MARS
  marsb <- earth(Y ~ ., data = na.omit(datab), degree = 2)
  
  # Extract alarms
  uni.alb <- mapply(">=", datab[,-1], 
    apply(marsb$cuts, 2, max))
  uni.alb[is.na(uni.alb)] <- TRUE
  alarmb <- apply(uni.alb, 1, all)
  alarms <- datab[alarmb,1]
  names(alarms) <- which(alarmb) 
  
  list(thresholds = apply(marsb$cuts, 2, max), alarms = alarms)
}

#' PRIM
PRIM.apply <- function(yb, xb){
  xb <- xb[complete.cases(yb),]
  yb <- na.omit(yb)
  n <- length(yb)
  
  # Peeling the box
  peelb <- peeling.sequence(yb, xb, alpha = .05, beta.stop = 6/n, 
    peeling.side = "left")
  
  # Extract box : box following the biggest increase in mean
  opt.sup <- which.max(diff(peelb$yfun) / -diff(peelb$support))
  supb <- peelb$support[opt.sup]
  chosenb <- extract.box(peelb, supb)
  boxb <- pasting.sequence(yb, xb, 
    small.box = chosenb$limits, peeling.side = "left", alpha = 0.01)
  
  # Extract alarms
  alarmb <- in.box(xb, boxb$limits)
  alarms <- yb[alarmb]
  names(alarms) <- which(alarmb)
    
  list(thresholds = boxb$limits[1,], alarms = alarms)
}

#' Classical
classical.apply <- function(yb, xb, extremes){
  xb <- xb[complete.cases(yb),]
  yb <- na.omit(yb)
  n <- length(yb)
  
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
  list(thresholds = unlist(thresholds), alarms = yb[alarmb]) 
}