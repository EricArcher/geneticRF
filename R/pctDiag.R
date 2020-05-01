#' @title Extract Percent Diagnosable
#' @description Calculate Random Forest confidence intervals
#' 
#' @param rf a \code{\link[randomForest]{randomForest}} object.
#' @param pd.vec a vector of of minimum vote values to return.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
pctDiag <- function(rf, pd.vec = NULL) {
  predicted <- apply(rf$votes, 1, function(x) colnames(rf$votes)[which.max(x)])
  vote.dfs <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    df <- data.frame(votes = votes, is.correct = predicted[y.i] == rf$y[y.i])
  }, simplify = F)
  min.p <-  1 / nlevels(rf$y)
  pd.vec <- if(is.null(pd.vec)) min.p else c(min.p, pd.vec)
  pd.vec <- pd.vec[pd.vec >= min.p]
  result <- sapply(pd.vec, function(p) {
    strata.pd <- sapply(vote.dfs, function(df) mean(df$is.correct & df$votes >= p))
    unname(strata.pd[which.min(strata.pd)])
  })
  names(result) <- paste("pd",  pd.vec, sep = ".")
  result
}
