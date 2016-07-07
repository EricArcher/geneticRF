#' @title Diagnosability Based on Shared Haplotypes
#' @description Maximum diagnosability based on shared haplotypes.
#' 
#' @param seq.df data.frame of stratified and aligned sequences from \code{\link{gtypes2rfDF}}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
sharedHapsPD <- function(seq.df) {
  hap.tbl <- table(makeHaps(seq.df[, -1, drop = F]), seq.df[, 1])
  min.freq.shared.haps <- rowSums(apply(hap.tbl, 1, function(y) {
    y[which.max(y)] <- 0
    y
  }))
  sum(min.freq.shared.haps) / sum(hap.tbl)
}
