#' @title Sequence Random Forest
#' @description Conduct Random Forest on stratified data.frame of sequences.
#' 
#' @param seq.df data.frame of stratified and aligned sequences from \code{\link{gtypes2rfDF}}.
#' @param replace sample with replacement? passed to \code{\link[randomForest]{randomForest}}.
#' @param sampsize number of samples to take from each strata. passed to \code{\link[randomForest]{randomForest}}.
#'   If \code{NULL}, value is set to the minumum sample size * \code{train.pct}.
#' @param train.pct proportion of strata to use for training if sampsize is \code{NULL}.
#' @param min.n minimum sample size to use. If sample size for any strata is less this, 
#'   this value will be used instead.
#' @param ... arguments passed to \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @import randomForest
#' @importFrom swfscMisc diversity
#' 
#' @export
#' 
sequenceRF <- function(seq.df, replace = FALSE, sampsize = NULL, train.pct = 0.5, min.n = 2, ...) {
  if(is.null(seq.df)) return(NULL)
  if(is.null(sampsize)) {
    strata.freq <- table(seq.df$strata)
    n <- ceiling(min(strata.freq) * train.pct)
    sampsize <- rep(n, length(strata.freq))
  }
  sampsize <- ifelse(sampsize < min.n, min.n, sampsize)
  rf <- randomForest(strata ~ ., data = seq.df, replace = replace, sampsize = sampsize, ...)
  
  overall.accuracy <- 1 - as.vector(rf$err.rate[nrow(rf$err.rate), "OOB"])
  conf.ci <- confIntRF(rf)
  min.diag <- which.min(conf.ci[-nrow(conf.ci), 1])
  diag.strata <- rownames(conf.ci)[min.diag]
  prior.diag <- as.vector(prop.table(table(rf$y))[diag.strata])
  hap.class.tbl <- classifyByHapFreq(seq.df)$class.tbl
  shared.hap.diag <- hap.class.tbl[diag.strata, "diagnosability"]
  haps <- makeHaps(seq.df[, -1, drop = FALSE])
  num.haps <- length(unique(haps))
  vs <- ncol(seq.df) - 1
  unique.hap.vs.vec <- uniqueHapBySite(seq.df)
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
    pd95 = pctDiag(rf, 0.95)[2],
    num.vs = vs,
    num.haps = num.haps,
    hap.div = swfscMisc::diversity(haps),
    num.vs.unique.hap = num.rf.vs.unique.hap,
    pct.vs.haps.shared = (vs - num.rf.vs.unique.hap) / vs,
    num.unique.haps.vs = num.unique.hap.rf.vs,
    pct.haps.vs.shared = (num.haps - num.unique.hap.rf.vs) / num.haps
  )
  rownames(smry) <- NULL
  
  list(smry = smry, rf = rf)
}
