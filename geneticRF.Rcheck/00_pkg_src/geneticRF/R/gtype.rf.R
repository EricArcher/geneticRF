#' @export gtype.rf
#' 
#' @title gtype Random Forest
#' @description Conduct Random Forest on a gtypes object.
#' 
#' @param g haploid \code{\link[strataG]{gtypes}} object with aligned sequences.
#' @param pairwise do analysis on all pairwise combinations of strata?
#' @param ... arguments passed to \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}). If \code{pairwise} is \code{TRUE} then 
#'   the \code{rf} element is a list of \code{randomForest} results for each row in
#'   \code{smry}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 

gtype.rf <- function(g, pairwise = FALSE, ...) {
  rf.func <- function(g) {
    seq.df <- create.seq.df(g)
    strata.freq <- table(seq.df$strata)
    min.n <- min(strata.freq)
    seq.df.rf(seq.df, replace = TRUE, sampsize = rep(min.n, length(strata.freq)), ...)
  }
  
  if(pairwise) {
    sp <- strata.pairs(g)
    result <- lapply(1:nrow(sp), function(i) {
      pair.g <- subset.gtypes(g, strata = unlist(sp[i, ]))
      gtype.rf(pair.g, ...)
    })
    list(smry = cbind(sp, do.call(rbind, lapply(result, function(x) x$smry))),
         rf = lapply(result, function(x) x$rf)
    )
  } else rf.func(g)
}
