context("test-mut_matrix")

# To test mut_matrix, we need to load the reference genome first.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# We re-use the data that is shipped with the package.
input <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

# Expected output
expected <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Run function
output <- mut_matrix(input, ref_genome)
output_longer <- mut_matrix(vcf_list = input, ref_genome = ref_genome, extension = 2)


# Perform tests

test_that("Output has correct class", {
  expect_true(inherits(output, "matrix"))
  expect_true(inherits(output_longer, "matrix"))
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output), c(96, 9))
  expect_equal(dim(output_longer), c(1536, 9))
})

test_that("Number of variants in output is correct", {
  expect_equal(colSums(output), elementNROWS(input))
  expect_equal(colSums(output_longer), elementNROWS(input))
})

test_that("transforms correctly", {
  expect_equal(output, expected)
})

test_that("a list is also acceptable input", {
  output_list <- mut_matrix(as.list(input), ref_genome)

  expect_equal(output_list, output)
  expect_equal(output_list, expected)
})

test_that("A single GR can also be used as input", {
  output_singlesample <- mut_matrix(input[[1]], ref_genome)
  expect_true(inherits(output_singlesample, "matrix"))
  expect_equal(dim(output_singlesample), c(96, 1))
})
