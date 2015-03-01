#' @export min.votes
#' 
#' @title Minimum Votes
#' @description Calculate Random Forest confidence intervals
#' 
#' @param rf a \code{\link[randomForest]{randomForest}} object.
#' @param mv.vec a vector of of minimum vote values to return.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 

min.votes <- function(rf, mv.vec) {
  vote.dfs <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    df <- data.frame(votes = votes, is.correct = rf$predicted[y.i] == rf$y[y.i])
    df[order(df$votes), ]
  }, simplify = F)
  
  sapply(1 - mv.vec, function(i) {
    min(sapply(vote.dfs, function(df) {
      cum.pct <- 1:nrow(df) / nrow(df)
      df$votes[min(which(cum.pct >= i))]
    }))
  })
}
