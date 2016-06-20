#' @title Make Haplotypes
#' @description Identify haplotypes in data.frame of aligned.
#' 
#' @param x data.frame of base pairs.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
makeHaps <- function(x) {
  dna.seq <- apply(x, 1, paste, collapse = "")
  most.freq <- sort(table(dna.seq), decreasing = T)
  as.numeric(factor(dna.seq, levels = names(most.freq)))
}
