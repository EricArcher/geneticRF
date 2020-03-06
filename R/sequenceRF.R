#' @title Sequence Random Forest
#' @description Conduct Random Forest on stratified data.frame of sequences.
#' 
#' @param x data.frame of stratified and aligned sequences from \code{\link{gtypes2rfDF}}.
#' @param replace sample with replacement? passed to \code{\link[randomForest]{randomForest}}.
#' @param sampsize number of samples to take from each strata. passed to \code{\link[randomForest]{randomForest}}.
#'   If \code{NULL}, value is set to the minumum sample size * \code{train.pct}.
#' @param train.pct proportion of strata to use for training if sampsize is \code{NULL}.
#' @param min.n minimum sample size to use. If sample size for any strata is less this, 
#'   this value will be used instead.
#' @param nrep number of rfPermute replicates.
#' @param conf.level confidence level for the \code{\link{binom.test}} confidence interval
#' @param ... arguments passed to \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list containing a data.frame of summary statistics (\code{smry}), and the 
#'   \code{randomForest} object (\code{rf}).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @seealso \code{\link[rfPermute]{classConfInt}}
#' 
#' @export
#' 
sequenceRF <- function(x, replace = FALSE, sampsize = NULL, 
                       train.pct = 0.5, min.n = 2, nrep = 0, 
                       conf.level = 0.95, ...) {
  if(is.null(x)) return(NULL)
  if(length(unique(x$stratum)) < 2) return(NULL)
  
  if(is.null(sampsize)) {
    strata.freq <- table(x$stratum)
    n <- ceiling(min(strata.freq) * train.pct)
    sampsize <- rep(n, length(strata.freq))
  }
  sampsize <- ifelse(sampsize < min.n, min.n, sampsize)
  rf <- rfPermute::rfPermute(
    stratum ~ ., data = x, replace = replace, 
    sampsize = sampsize, nrep = nrep, ...
  )
  
  overall.accuracy <- 1 - as.vector(rf$err.rate[nrow(rf$err.rate), "OOB"])
  ci <- rfPermute::confusionMatrix(rf, conf.level = conf.level)
  ci <- ci[, (length(rf$classes) + 1):(length(rf$classes) + 3)]
  min.diag <- which.min(ci[-nrow(ci), 1])
  diag.strata <- rownames(ci)[min.diag]
  hap.class.tbl <- classifyByHapFreq(x)$class.tbl
  shared.hap.diag <- hap.class.tbl[diag.strata, "diagnosability"]
  haps <- makeHaps(x[, -1, drop = FALSE])
  num.haps <- length(unique(haps))
  vs <- ncol(x) - 1
  unique.hap.vs.vec <- uniqueHapBySite(x)
  num.rf.vs.unique.hap <- length(unique.hap.vs.vec)
  num.unique.hap.rf.vs <- length(table(unique.hap.vs.vec))
  
  smry <- data.frame(
    overall.accuracy = overall.accuracy * 100,
    diag.strata = diag.strata,
    diagnosability = ci[min.diag, 1],
    diag.lci = ci[min.diag, 2],
    diag.uci = ci[min.diag, 3],
    shared.hap.diag = shared.hap.diag,
    pd95 = pctDiag(rf, 0.95)[2],
    num.vs = vs,
    num.haps = num.haps,
    hap.div = swfscMisc::diversity(haps)
  )
  rownames(smry) <- NULL
  
  list(smry = smry, rf = rf)
}
