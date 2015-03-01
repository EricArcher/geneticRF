#' @export create.seq.df
#' @import strataG
#' 
#' @title Create Sequence Data.Frame
#' @description Create data.frame of variable sites from gtypes object.
#' 
#' @param g haploid \code{\link[strataG]{gtypes}} object with aligned sequences.
#' @param label label to add to beginning of each stratum name.
#' 
#' @return a \code{data.frame} where the first column lists the (\code{strata}) and every column
#'   column afterwards is a variable site. All columns are factors.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 

create.seq.df <- function(g, label = NULL) {
  stopifnot.gtypes(g, "haploid")
  stopifnot.aligned(g$sequences)
  options(stringsAsFactors = F)

  # extract strata data frame
  g <- decode(g)
  strata.df <- as.data.frame(g$genotypes)
  strata.df$ids <- rownames(strata.df)

  # remove rows without stratifications
  strata.df <- na.omit(strata.df)

  # sequences to use are those in the updated strata data.frame
  dna.seqs <- g$sequences[which(names(g$sequences) %in% strata.df$haplotype)]

  # extract variable sites for these sequences and create sequence matrix
  var.sites <- variable.sites(dna.seqs)
  var.seq.mat <- tolower(do.call(rbind, var.sites$sites))
  sites <- paste("site", colnames(var.sites$site.freqs), sep = ".")

  # only use sites with valid bases
  to.keep <- apply(var.seq.mat, 2, function(x) all(x %in% c("a", "c", "g", "t", "-", ".")))
  if(sum(to.keep) == 0) return(NULL)
  var.seq.mat <- cbind(var.seq.mat[, to.keep])
  colnames(var.seq.mat) <- sites[to.keep]

  # create factors of variable site columns
  var.seq.df <- do.call(data.frame, lapply(colnames(var.seq.mat), function(x) factor(var.seq.mat[, x])))
  colnames(var.seq.df) <- colnames(var.seq.mat)
  var.seq.df <- cbind(haplotype = names(var.sites$sites), var.seq.df)
  seq.df <- merge(strata.df, var.seq.df, all = T, by = "haplotype")

  # remove any rows with missing data and remove id and haplotype columns
  seq.df <- na.omit(seq.df)
  seq.df$strata <- factor(seq.df$strata)
  rownames(seq.df) <- seq.df$ids
  seq.df$id <- seq.df$haplotype <- NULL
  sites <- seq.df[, -1, drop = FALSE]

  # remove sites where substitutions are only represented by one individual
  to.keep <- apply(sites, 2, function(x) sum(table(x) > 1) > 1)
  if(sum(to.keep) == 0) return(NULL)

  strata <- if(is.null(label)) seq.df$strata else paste(label, seq.df$strata)
  cbind(strata = factor(strata), sites[, to.keep, drop = FALSE])
}
