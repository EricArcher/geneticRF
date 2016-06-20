#' @title Create Sequence Data.Frame
#' @description Create data.frame of variable sites from gtypes object.
#' 
#' @param g haploid \code{\link[strataG]{gtypes}} object with aligned sequences.
#' @param gene number or name of gene to use from multidna \code{@sequences} slot.
#' @param label label to add to beginning of each stratum name.
#' 
#' @return a \code{data.frame} where the first column lists the (\code{strata}) and every column
#'   column afterwards is a variable site. All columns are factors.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @import strataG
#' @export

gtypes2seqDF <- function(g, gene = 1, label = NULL) {
  if(ploidy(g) != 1) stop("'g' should be a haploid gtypes object.")
  
  # extract stratified sequences
  df <- na.omit(as.data.frame(g))
  colnames(df)[3] <- "haplotype"
  dna.seqs <- sequences(g[df$id, , , drop = TRUE], gene)
  
  # extract variable sites for these sequences and create sequence matrix
  var.sites <- variableSites(dna.seqs)
  var.seq.mat <- tolower(do.call(rbind, as.character(as.list(var.sites$sites))))
  sites <- paste("site", colnames(var.sites$site.freqs), sep = ".")
  
  # only use sites with valid bases
  to.keep <- apply(var.seq.mat, 2, function(x) all(x %in% c("a", "c", "g", "t", "-", ".")))
  if(sum(to.keep) == 0) return(NULL)
  var.seq.mat <- cbind(var.seq.mat[, to.keep])
  colnames(var.seq.mat) <- sites[to.keep]
  
  # create factors of variable site columns
  seq.df <- do.call(data.frame, lapply(colnames(var.seq.mat), function(x) factor(var.seq.mat[, x])))
  colnames(seq.df) <- colnames(var.seq.mat)
  
  # add strata and remove any rows with missing data
  seq.df <- data.frame(cbind(strata = df$strata, seq.df[df$haplotype, ]))
  rownames(seq.df) <- rownames(df)
  seq.df <- na.omit(seq.df)
  
  # remove sites where substitutions are only represented by one individual
  sites <- seq.df[, -1, drop = FALSE]
  to.keep <- apply(sites, 2, function(x) sum(table(x) > 1) > 1)
  if(sum(to.keep) == 0) return(NULL)
  
  st <- if(is.null(label)) seq.df$strata else paste(label, seq.df$strata)
  cbind(strata = factor(st), sites[, to.keep, drop = FALSE])
}
