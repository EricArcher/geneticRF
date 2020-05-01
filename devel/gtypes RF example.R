rm(list = ls())
library(strataG)
library(geneticRF)
library(rfPermute)
library(ape)
data(dolph.strata)
data(dolph.haps)

# Write example data ------------------------------------------------------

write.csv(dolph.strata[, c(1, 4, 2)], file = "strata.csv", row.names = FALSE)
dloop.seqs <- as.DNAbin(dolph.haps)
write.fasta(dloop.seqs, file = "dloop.fasta")
rm(list = ls())

# Read data ---------------------------------------------------------------

ex.strata <- read.csv("strata.csv", stringsAsFactors = FALSE)
ex.seqs <- read.fasta("dloop.fasta")

# Create gtypes -----------------------------------------------------------

ex.g <- df2gtypes(ex.strata, ploidy = 1, sequences = ex.seqs)
summary(ex.g)

# Run Random Forest -------------------------------------------------------

ex.rf <- gtypesRF(ex.g)
confusionMatrix(ex.rf$rf)
