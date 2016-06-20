rm(list = ls())
library(strataG)
library(geneticRF)

data(dolph.strata)
data(dolph.haps)
mtdna <- gtypes(dolph.strata, id.col = 1, strata.col = 3, locus.col = 4, dna.seq = dolph.haps)
summary(mtdna)

gtype.rf(mtdna)

rf.species.id(mtdna, c("Coastal", "Offshore.North"), "Offshore.South")
