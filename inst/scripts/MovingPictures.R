# -----------------------------------------------------------------------------
# MovingPictures.R
#
# Preparation script for the miaTime example dataset `MovingPictures`
# (see `?miaTime::MovingPictures`).
#
# Source: the phyloseq object `MovingPictures` of the MicrobeDS package
# (https://github.com/twbattaglia/MicrobeDS, CC0), from Qiita study 550
# (Caporaso et al. 2011, Genome Biology 12:R50).
#
# The OTUs are agglomerated to genus to keep the package small. OTUs without
# a genus are kept at their lowest assigned rank, so no reads are lost.
#
# Run this script from the package root.
# -----------------------------------------------------------------------------

library(mia)

file <- tempfile(fileext = ".rda")
download.file(
    "https://raw.githubusercontent.com/twbattaglia/MicrobeDS/master/data/MovingPictures.rda",
    file, mode = "wb")
load(file)

tse <- convertFromPhyloseq(MovingPictures)
tse <- agglomerateByRank(
    tse, rank = "Genus", empty.rm = FALSE, update.tree = TRUE)

MovingPictures <- tse
save(MovingPictures, file = "data/MovingPictures.rda", compress = "xz")
