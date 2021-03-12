#########################################################
#            Main code to produce results in
#                                                  
#     "Objective methods for threshold assessment in 
#          weather-health watch warning systems"
#
#             Author: Pierre Masselot
#
#                    Journal
#
#              Bootstrap simulations
#
#########################################################

#--------------------------------
# Perform resampling
#--------------------------------

# Number of Bootstrap samples
B <- 1000

# Create data blocks
datablock <- split(1:n, datatab$Year)
bpool <- 1:length(datablock)

# Draw bootstrap samples
set.seed(12345)
bsamples <- replicate(B, sample(datablock, replace = T), simplify = F)
bsamples <- lapply(bsamples, unlist)

#--------------------------------
# Prepare parallel computing
#--------------------------------

cl <- makeCluster(detectCores() - 2)
# Transfer objects in cluster
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(quantmod)
  library(earth)
  library(splines)
  library(lattice)
  library(primr)
  library(mgcv)
  library(partykit)
  library(AIM)
})

#--------------------------------
# Apply algorithms
#--------------------------------

# Prepare object containing results
bootRes <- list()

#----- MOB -----
bootRes[["MOB"]] <- parLapply(cl, bsamples, function(b){
  dat <- datatab[b,]
  MOB.apply(yb = dat$Death, xb = dat[,indicator.names], 
    zb = cbind(dat$Tmoy, ns(dat$dos, 4), ns(dat$Year, round(nyear / 10))),
    family = "quasipoisson", minsize = 10)
}) 

#----- MARS -----
bootRes[["MARS"]] <- parLapply(cl, bsamples, function(b){
  dat <- datatab[b,]
  MARS.apply(yb = dat$Death, xb = dat[,indicator.names], 
    zb = dat[,c("dos", "Year")], degree = 2, endspan = 10, 
    glm = list(family = "quasipoisson"))
})

#----- PRIM -----
bootRes[["PRIM"]] <- parLapply(cl, bsamples, function(b){
  dat <- datatab[b,]
  PRIM.apply(yb = dat$Death, xb = dat[,indicator.names], 
    zb = cbind(dat$Tmoy, ns(dat$dos, 4), ns(dat$Year, round(nyear / 10))),
    RRind = 1, beta.stop = 10/n, peeling.side = -1, family = "quasipoisson")    
})

#----- AIM -----
bootRes[["AIM"]] <- parLapply(cl, bsamples, function(b){
  dat <- datatab[b,]
  AIM.apply(yb = dat$Death, xb = dat[,indicator.names], 
    zb = dat[,c("dos", "Year")], backfit = T, mincut = 10/n)
})

# Stop parallel
stopCluster(cl)

# Save results
save(bsamples, bootRes, file = "Results/21_BootRes.RData")