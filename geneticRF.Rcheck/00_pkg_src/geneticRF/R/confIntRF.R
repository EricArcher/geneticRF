#' @title Confidence Intervals
#' @description Calculate Random Forest confidence intervals
#' 
#' @param rf a \code{\link[randomForest]{randomForest}} object
#' @param conf.level confidence level for the \code{\link{binom.test}} confidence interval
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
confIntRF <- function(rf, conf.level = 0.95) {
  conf <- rf$confusion
  result <- t(apply(conf, 1, function(y) {
    pct.correct <- as.vector(1 - y[length(y)])
    total.n <- sum(y[-length(y)])
    ci <- binom.test(round(total.n * pct.correct, 0), total.n, pct.correct, conf.level = conf.level)$conf.int
    names(ci) <- paste(c("LCI", "UCI"), conf.level, sep = "_")
    c(pct.correct = pct.correct, ci)
  }))
  total.n <- sum(conf[, -ncol(conf)])
  num.correct <- sum(diag(conf[, -ncol(conf)]))
  pct.correct <- num.correct / total.n
  ci <- binom.test(num.correct, total.n, pct.correct, conf.level = conf.level)$conf.int
  rbind(result, Total = c(pct.correct, ci))
}
