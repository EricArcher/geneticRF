#' @title gtype Random Forest
#' @description Conduct Random Forest on a gtypes object.
#' 
#' @param g haploid \code{\link[strataG]{gtypes}} object with aligned sequences.
#' @param gene number or name of gene to use from multidna \code{@sequences} slot.
#' @param pairwise do analysis on all pairwise combinations of strata?
#' @param ... arguments passed to \code{\link{sequenceRF}} and \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}). If \code{pairwise} is \code{TRUE} then 
#'   the \code{rf} element is a list of \code{randomForest} results for each row in
#'   \code{smry}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export

gtypesRF <- function(g, gene = 1, pairwise = FALSE, ...) {
  # check for and set names to sampsize
  args <- list(...)
  ss.match <- which(pmatch(names(args), "sampsize") == 1)
  if(length(ss.match) == 1) {
    temp.ss <- args[[ss.match]]
    if(is.null(names(temp.ss))) {
      names(temp.ss) <- levels(strata(g))
      args[[ss.match]] <- temp.ss
      names(args)[ss.match] <- "sampsize"
    }
  }
  
  if(pairwise) {
    sp <- strataG:::.strataPairs(g)
    result <- lapply(1:nrow(sp), function(i) {
      pair.g <- g[, , unlist(sp[i, ])]
      do.call(gtypesRF, c(list(g = pair.g, pairwise = FALSE), args))
      #gtypesRF(pair.g, pairwise = FALSE, ...)
    })
    list(smry = cbind(sp, do.call(rbind, lapply(result, function(x) x$smry))),
         rf = lapply(result, function(x) x$rf)
    )
  } else {
    args <- c(list(seq.df = gtypes2seqDF(g, gene = gene)), args)
    do.call(sequenceRF, args)
    #sequenceRF(seq.df, ...)
  }
}
