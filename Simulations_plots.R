#########################################################
#
#                  Simulation plots
#                                                  
#########################################################

#---------------------------------------------------
#                 Parameters
#---------------------------------------------------

result.path <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 1 - thresholds/Simulations"
setwd(result.path)

comparison_name <- "Basic_case"

simus <- c("Basic1", "Basic3", "Basic5", "Basic10")
pars <- c(1, 3, 5, 10)

models <- c("CART", "MARS", "PRIM")
mod.names <- c("CART", "MARS", "PRIM")

sref <- 

#--- Deriving parameters ---
ns <- length(simus)
nm <- length(models)

#---------------------------------------------------
#                 Reading results
#---------------------------------------------------

result_list <- vector("list", nm)
names(result_list) <- mod.names

for (s in 1:ns){
  load(sprintf("%s/All_results.RData", simus[s]))
  results <- results[models]
  
  ##### Thresholds
  thresholds <- lapply(results, sapply, "[[", "thresholds")
  
  # Bias computing 
  meanThresh <- sapply(thresholds, apply, 1, mean, na.rm = T)
  bias <- meanThresh - 
    matrix(s, nrow = nrow(meanThresh), ncol = length(results))
}