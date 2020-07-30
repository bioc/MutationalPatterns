context("test-bin_mutation_density")


# Read grl
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

## Load the corresponding reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

## Determine region density
output <- bin_mutation_density(grl, ref_genome, nrbins = 3)


# Use manual cutoffs
output_man <- bin_mutation_density(grl, ref_genome, man_dens_cutoffs = c(0, 2e-08, 1))


# Tests
test_that("Output has correct class", {
  expect_true(inherits(output, "CompressedGRangesList"))
  expect_true(inherits(output_man, "CompressedGRangesList"))
})

test_that("Output has correct dimensions", {
  expect_equal(length(output), 3)
  expect_equal(length(output_man), 2)
  expect_equal(as.vector(S4Vectors::elementNROWS(output)), c(30, 11, 2))
  expect_equal(as.vector(S4Vectors::elementNROWS(output_man)), c(25, 4))
})
