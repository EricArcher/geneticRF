#' @title Make Haplotypes
#' @description Identify haplotypes in data.frame of aligned.
#' 
#' @param x data.frame of base pairs.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
uniqueHapBySite <- function(x) {
  haps <- makeHaps(x)
  result <- unlist(apply(x[, -1, drop = F], 2, function(b) {
    site.tbl <- table(haps, b)
    shared.site <- apply(site.tbl, 2, function(y) sum(y > 0))
    j <- which(shared.site == 1)
    if(length(j) > 0) {
      rownames(site.tbl)[which(site.tbl[, j] > 0)]
    } else NA
  }))
  result[!is.na(result)]
}
