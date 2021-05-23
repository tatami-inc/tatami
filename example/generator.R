# This writes a hypothetical SummarizedExperiment's contents into a several binary blobs.

dir.create("colData")
writeBin(1:26, "colData/stuff", endian="big")

dir.create("assays")
dir.create('assays/counts')
writeBin(c(10L, 26L), "assays/counts/dim", endian="big")
writeBin(500 + 1:260, "assays/counts/values", endian="big")
