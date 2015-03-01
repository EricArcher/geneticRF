#' @export seq.df.rf
#' @import randomForest
#' 
#' @title Sequence Random Forest
#' @description Conduct Random Forest on stratified data.frame of sequences.
#' 
#' @param seq.df data.frame of stratified and aligned sequences from \code{\link{create.seq.df}}.
#' @param ... arguments passed to \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 

seq.df.rf <- function(seq.df, ...) {
  if(is.null(seq.df)) return(NULL)
  rf <- randomForest(strata ~ ., data = seq.df, ...)
  
  overall.accuracy <- 1 - as.vector(rf$err.rate[nrow(rf$err.rate), "OOB"])
  conf.ci <- rf.conf.int(rf)
  min.diag <- which.min(conf.ci[, 1])
  diag.strata <- rownames(conf.ci)[min.diag]
  prior.diag <- as.vector(prop.table(table(rf$y))[diag.strata])
  hap.class.tbl <- classify.by.hap.freq(seq.df)$class.tbl
  shared.hap.diag <- hap.class.tbl[diag.strata, "diagnosability"]
  haps <- make.haps(seq.df[, -1, drop = FALSE])
  num.haps <- length(unique(haps))
  vs <- ncol(seq.df) - 1
  unique.hap.vs.vec <- unique.hap.by.site(seq.df)
  num.rf.vs.unique.hap <- length(unique.hap.vs.vec)
  num.unique.hap.rf.vs <- length(table(unique.hap.vs.vec))

  smry <- data.frame(
    overall.accuracy = overall.accuracy,
    diag.strata = diag.strata,
    diagnosability = conf.ci[min.diag, 1],
    diag.lci = conf.ci[min.diag, 2],
    diag.uci = conf.ci[min.diag, 3],
    prior.diag = prior.diag,
    shared.hap.diag = shared.hap.diag,
    pd95 = pct.diag(rf, 0.95)[2],
    num.vs = vs,
    num.haps = num.haps,
    hap.div = diversity(haps),
    num.vs.unique.hap = num.rf.vs.unique.hap,
    pct.vs.haps.shared = (vs - num.rf.vs.unique.hap) / vs,
    num.unique.haps.vs = num.unique.hap.rf.vs,
    pct.haps.vs.shared = (num.haps - num.unique.hap.rf.vs) / num.haps
  )
  rownames(smry) <- NULL
  
  list(smry = smry, rf = rf)
}
