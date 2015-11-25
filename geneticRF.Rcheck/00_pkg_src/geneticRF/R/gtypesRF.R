#' @title gtype Random Forest
#' @description Conduct Random Forest on a gtypes object.
#' 
#' @param g haploid \code{\link[strataG]{gtypes}} object with aligned sequences.
#' @param gene number or name of gene to use from multidna \code{@sequences} slot.
#' @param pairwise do analysis on all pairwise combinations of strata?
#' @param replace sample with replacement? passed to \code{\link[randomForest]{randomForest}}.
#' @param sampsize number of samples to take from each strata. passed to \code{\link[randomForest]{randomForest}}.
#' @param ... arguments passed to \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}). If \code{pairwise} is \code{TRUE} then 
#'   the \code{rf} element is a list of \code{randomForest} results for each row in
#'   \code{smry}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export

gtypesRF <- function(g, gene = 1, pairwise = FALSE, replace = FALSE, sampsize = NULL, ...) {
  rf.func <- function(g, replace = FALSE, sampsize = NULL) {
    seq.df <- gtypes2seqDF(g, gene = gene)
    strata.freq <- table(seq.df$strata)
    min.n <- min(strata.freq)
    if(is.null(sampsize)) sampsize <- rep(min.n, length(strata.freq))
    sequenceRF(seq.df, replace = replace, sampsize = sampsize, ...)
  }
  
  if(!is.null(sampsize)) names(sampsize) <- levels(strata(g))
  if(pairwise) {
    sp <- strataG:::.strataPairs(g)
    result <- lapply(1:nrow(sp), function(i) {
      pair.g <- g[, , unlist(sp[i, ])]
      gtypesRF(pair.g, replace = replace, sampsize = sampsize[levels(strata(pair.g))], ...)
    })
    list(smry = cbind(sp, do.call(rbind, lapply(result, function(x) x$smry))),
         rf = lapply(result, function(x) x$rf)
    )
  } else rf.func(g, replace = replace, sampsize = sampsize)
}
