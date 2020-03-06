#' @title Random Forest Species ID
#' @description Predict species for unknown samples based on reference samples using Random Forest.
#' 
#' @param g haploid \code{\link[strataG]{gtypes}} object with aligned sequences.
#' @param ref.strata a character vector of 2 or more strata in \code{g} top use as reference samples.
#' @param unk.strata a character vector of strata representing unknown samples.
#' @param ... arguments passed to \code{\link{sequenceRF}} and \code{\link[randomForest]{randomForest}}.
#' 
#' @return a list with the following elements:
#' \tabular{ll}{
#'   \code{pred} \tab vector of species predictions.\cr
#'   \code{prob} \tab matrix of species assignment probabilities.\cr
#'   \code{rf} \tab the \code{randomForest} model object.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @export
#' 
rfSpeciesID <- function(g, ref.strata = NULL, unk.strata = NULL, ...) {
  g.strata <- sort(strataG::getStrataNames(g))
  if(is.null(ref.strata)) ref.strata <- setdiff(g.strata, unk.strata)
  if(is.null(unk.strata)) unk.strata <- setdiff(g.strata, ref.strata)
  if(length(intersect(ref.strata, unk.strata)) > 0) stop("'ref.strata' and 'unk.strata' cannot overlap.")
  if(length(unique(ref.strata)) < 2) stop("'ref.strata' must have 2 or more strata.")
  sub.g <- g[, , c(ref.strata, unk.strata), drop = TRUE]
  seq.df <- gtypes2rfDF(sub.g)
  ref.df <- seq.df[seq.df$stratum %in% ref.strata, ]
  ref.df$stratum <- droplevels(ref.df$stratum)
  ref.rf <- sequenceRF(ref.df, ...)
  if(length(unk.strata) == 0) return(list(pred = NA, prob = NA, rf = ref.rf))
  unk.df <- seq.df[seq.df$stratum %in% unk.strata, ]
  unk.df$stratum <- droplevels(unk.df$stratum)
  unk.pred <- stats::predict(ref.rf$rf, unk.df, type = "response")
  unk.prob <- stats::predict(ref.rf$rf, unk.df, type = "prob")
  list(pred = unk.pred, prob = unk.prob, rf = ref.rf)
}
