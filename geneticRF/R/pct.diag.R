#' @export pct.diag
#' 
#' @title Extract Percent Diagnosable
#' @description Calculate Random Forest confidence intervals
#' 
#' @param rf a \code{\link[randomForest]{randomForest}} object.
#' @param pd.vec a vector of of minimum vote values to return.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 

pct.diag <- function(rf, pd.vec) {
  vote.dfs <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    df <- data.frame(votes = votes, is.correct = rf$predicted[y.i] == rf$y[y.i])
    df[order(df$votes), ]
  }, simplify = F)
  strata.pd <- sapply(vote.dfs, function(df) sum(df$is.correct) / nrow(df))
  is.min <- which.min(strata.pd)
  min.df <- vote.dfs[[is.min]]
  result <- sapply(pd.vec, function(i) sum(min.df$votes > i) / nrow(min.df))
  result <- c(strata.pd[is.min], result)
  names(result) <- paste("pd", c("0.5", pd.vec), sep = ".")
  result
}
