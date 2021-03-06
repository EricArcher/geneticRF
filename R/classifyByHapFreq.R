#' @title Classify By Haplotype Frequency
#' @description Classify samples by haplotype frequency.
#' 
#' @param seq.df data.frame of stratified and aligned sequences from \code{\link{gtypes2rfDF}}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
classifyByHapFreq <- function(seq.df) {
  haps <- makeHaps(seq.df[, -1, drop = F])
  hap.freq <- prop.table(table(haps, seq.df$stratum), 1)
  pred <- sapply(haps, function(i) colnames(hap.freq)[which.max(hap.freq[as.character(i), ])])
  strata <- unique(seq.df$stratum)
  class.tbl <- table(factor(seq.df$stratum, levels = strata), factor(pred, levels = strata))
  overall.diag <- sum(diag(class.tbl)) / sum(class.tbl)
  class.tbl <- cbind(class.tbl, diagnosability = diag(class.tbl) / rowSums(class.tbl))
  list(pred = pred, class.tbl = class.tbl, overall.diag = overall.diag)
}
