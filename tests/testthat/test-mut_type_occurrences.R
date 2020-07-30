context("mut_type_occurrences")

# Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

## Load a reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Get the type occurrences for all VCF objects.
output <- mut_type_occurrences(vcfs, ref_genome)

# Get type occurence for single sample
output_single_sample <- mut_type_occurrences(vcfs[[1]], ref_genome)

# Get type occurence for few muts
output_fewmuts <- mut_type_occurrences(vcfs[[1]][1:2], ref_genome)

test_that("Output has correct class", {
  expect_true(inherits(output, "data.frame"))
  expect_true(inherits(output_single_sample, "data.frame"))
  expect_true(inherits(output_fewmuts, "data.frame"))
})

test_that("Outpus has correct dimensions", {
  expect_equal(dim(output), c(9, 8))
  expect_equal(dim(output_single_sample), c(1, 8))
  expect_equal(dim(output_fewmuts), c(1, 8))
})

test_that("Transforms correctly", {
  expect_equal(
    output_single_sample,
    structure(list(
      `C>A` = 28L, `C>G` = 5L, `C>T` = 109L, `T>A` = 12L,
      `T>C` = 30L, `T>G` = 12L, `C>T at CpG` = 59L, `C>T other` = 50L
    ),
    row.names = "My_sample", class = "data.frame"
    )
  )
})
