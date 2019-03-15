#' F-score
#'
#' Computes the F-score using observed and predicted indices of the event.
#'
#' @param predicted Vector of index of predicted events.
#' @param observed Vector of index of observed events.
#' @param beta Beta value of the F-score. If larger than 1, more importance
#'  is given to recall while if lower, more importance is given to precision.
#'
#' @detail The F-score is based on recall which is identical to sensivity,
#'  i.e. the proportion of events predicted, and sensitivity which is the 
#'  proportion of predicted events that are true events. 
#'
#' @references
#'
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