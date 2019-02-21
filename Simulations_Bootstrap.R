#########################################################
#            Bootstrap simulations in paper
#                                                  
#             follow-up to the main code
#########################################################

#---------------------------------------------------
#                 Initialization
#---------------------------------------------------
# parameters
B <- 100

# Object containing results
boot.thresh <- boot.alarms <- list()
cols <- c()

# Formula
rhs <- paste(sprintf("Z%s", which.X), collapse = " + ")
bform <- as.formula(sprintf("OM ~ %s", rhs))

#---------------------------------------------------
#            Create bootstrap samples
#---------------------------------------------------

# Split data in blocks
datablock <- split(datatab, datatab$Date$year)
bpool <- 1:length(datablock)

# Draw bootstrap samples
binds <- replicate(B, sample(bpool, replace = T))

#---------------------------------------------------
#                 CART application
#---------------------------------------------------

cols["CART"] <- "forestgreen"

# Object containing results
boot.thresh[["CART"]] <- matrix(NA, nrow = B, ncol = p, 
  dimnames = list(NULL, which.X))
boot.alarms[["CART"]] <- vector("list", B)

# Loop
for(b in 1:B){
  # Construct the bth sample
  datab <- datablock[binds[,b]]
  datab <- Reduce(rbind, datab)

  # Grow the tree
  treeb <- rpart(bform, data = datab, method = "anova", 
    control = rpart.control(minsplit = 10, cp = 0.0001))
  minind <- which.min(treeb$cptable[,"xerror"])
  best_cp <- treeb$cptable[minind, "CP"]
  treeb <- prune(treeb, best_cp)
  
  # Extract thresholds
  boxb <- get.box(treeb)
  boot.thresh[["CART"]][b,] <- boxb$box[1,]
  
  # Alarms
  uni.alb <- mapply(">=", datab[,sprintf("Z%s", which.X)], 
    boot.thresh[["CART"]][b,])
  alarmb <- apply(uni.alb, 1, all) 
  boot.alarms[["CART"]][[b]] <- datab$OM[alarmb]
}

#---------------------------------------------------
#                 MARS application
#---------------------------------------------------

cols["MARS"] <- "cornflowerblue"

# Object containing results
boot.thresh[["MARS"]] <- matrix(NA, nrow = B, ncol = p, 
  dimnames = list(NULL, which.X))
boot.alarms[["MARS"]] <- vector("list", B)
  
# Loop
for(b in 1:B){
  # Construct the bth sample
  datab <- datablock[binds[,b]]
  datab <- Reduce(rbind, datab)
  
  # Apply MARS
  marsb <- earth(bform, data = na.omit(datab))
  
  # Extract thresholds
  boot.thresh[["MARS"]][b,] <- apply(marsb$cuts, 2, max)
  
  # Extract alarms
  uni.alb <- mapply(">=", datab[,sprintf("Z%s", which.X)], 
    boot.thresh[["MARS"]][b,])
  alarmb <- apply(uni.alb, 1, all) 
  boot.alarms[["MARS"]][[b]] <- datab$OM[alarmb]
}

#---------------------------------------------------
#                 PRIM application
#---------------------------------------------------

cols["PRIM"] <- "firebrick"

# Object containing results
boot.thresh[["PRIM"]] <- matrix(NA, nrow = B, ncol = p, 
  dimnames = list(NULL, which.X))
boot.alarms[["PRIM"]] <- vector("list", B)
  
# Loop
for(b in 1:B){
  # Construct the bth sample
  datab <- datablock[binds[,b]]
  datab <- Reduce(rbind, datab)
  
  # Peeling the box
  peelb <- peeling.sequence(datab[,"OM"], datab[,sprintf("Z%s", which.X)], 
    alpha = .05, beta.stop = 6/n, peeling.side = "left")
  
  # Extract box : box following the biggest increase in mean
  supb <- peelb$support[which.max(diff(peelb$yfun))]
  chosenb <- extract.box(peelb, supb)
  boxb <- pasting.sequence(datab[,"OM"], datab[,sprintf("Z%s", which.X)], 
    small.box = chosenb$limits, peeling.side = "left", alpha = 0.01)
  
  # Extract thresholds
  boot.thresh[["PRIM"]][b,] <- boxb$limits[1,]
  
  # Extract alarms
  boot.alarms[["PRIM"]][[b]] <- in.box(datab[,sprintf("Z%s", which.X)], 
    boxb$limits, datab$OM)
}

#---------------------------------------------------
#                     Plots
#---------------------------------------------------  

# Thresholds by indicator
indic.thresholds <- list()
for (j in 1:p){
  indic.thresholds[[sprintf("Z%s", which.X[j])]] <- sapply(boot.thresh, "[", ,j)
}

# Plot thresholds
x11(title = "Thresholds", height = 5, width = 10)
par(mfrow = c(1, p))
for (j in 1:p){
  boxplot(indic.thresholds[[j]], border = cols, lwd = 2, ylab = "Threshold", 
    main = sprintf("Z%s", which.X[j]), cex.lab = 1.3, cex.axis = 1.2, xlab = "",
    varwidth = T)
}

# Number of alarms
nalarms <- sapply(boot.alarms, sapply, length)
x11(title = "Alarm number")
boxplot(nalarms, border = cols, lwd = 2, ylab = "Threshold", xlab = "", 
  cex.lab = 1.3, cex.axis = 1.2, varwidth = T) 
  
# Mean value of alarms
nalarms <- sapply(boot.alarms, sapply, mean)
x11(title = "Alarm mean")
boxplot(nalarms, border = cols, lwd = 2, ylab = "Threshold", xlab = "", 
  cex.lab = 1.3, cex.axis = 1.2, varwidth = T) 