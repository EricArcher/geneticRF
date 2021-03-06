% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtypes2rfDF.R
\name{gtypes2rfDF}
\alias{gtypes2rfDF}
\title{Create Sequence Data.Frame}
\usage{
gtypes2rfDF(g, gene = 1, label = NULL)
}
\arguments{
\item{g}{a \code{\link[strataG]{gtypes}} object. If haploid, it must have 
aligned sequences.}

\item{gene}{number or name of gene to use from multidna \code{@sequences} slot.}

\item{label}{label to add to beginning of each stratum name.}
}
\value{
a \code{data.frame} where the first column lists the (\code{strata}) and every column
  column afterwards is a variable site. All columns are factors.
}
\description{
Create data.frame of variable sites from gtypes object.
}
\author{
Eric Archer \email{eric.archer@noaa.gov}
}
