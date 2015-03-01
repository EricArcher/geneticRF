#' @export collapse.to.haps
#' 
#' @title Collapse To Haplotypes
#' @description Collapse sequence data.frame to unique haplotypes.
#' 
#' @param seq.df data.frame of stratified and aligned sequences from \code{\link{create.seq.df}}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 

collapse.to.haps <- function(seq.df) {
  haps <- make.haps(seq.df[, -1, drop = F])
  seq.df[!duplicated(haps), , drop = F]
}
