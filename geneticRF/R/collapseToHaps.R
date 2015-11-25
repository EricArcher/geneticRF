#' @title Collapse To Haplotypes
#' @description Collapse sequence data.frame to unique haplotypes.
#' 
#' @param seq.df data.frame of stratified and aligned sequences from \code{\link{gtypes2seqDF}}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
collapseToHaps <- function(seq.df) {
  haps <- makeHaps(seq.df[, -1, drop = F])
  seq.df[!duplicated(haps), , drop = F]
}
