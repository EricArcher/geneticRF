#' @title Classify By Haplotype Frequency
#' @description Classify samples by haplotype frequency.
#' 
#' @param seq.df data.frame of stratified and aligned sequences from \code{\link{gtypes2seqDF}}.
#' @param ... arguments passed to \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
classifyByHapFreq <- function(seq.df) {
  haps <- makeHaps(seq.df[, -1, drop = F])
  hap.freq <- prop.table(table(haps, seq.df$strata), 1)
  pred <- sapply(haps, function(i) colnames(hap.freq)[which.max(hap.freq[as.character(i), ])])
  strata <- unique(seq.df$strata)
  class.tbl <- table(factor(seq.df$strata, levels = strata), factor(pred, levels = strata))
  overall.diag <- sum(diag(class.tbl)) / sum(class.tbl)
  class.tbl <- cbind(class.tbl, diagnosability = diag(class.tbl) / rowSums(class.tbl))
  list(pred = pred, class.tbl = class.tbl, overall.diag = overall.diag)
}
