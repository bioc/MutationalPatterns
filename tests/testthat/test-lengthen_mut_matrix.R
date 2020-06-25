context("test-lengthen_mut_matrix")

# Read in mut_matrix
input <- readRDS(system.file("states/mut_mat_splitregions.rds",
  package = "MutationalPatterns"
))

# Read in indel
input_indel <- readRDS(system.file("states/blood_indels_counts_split_region.rds",
  package = "MutationalPatterns"
))

## Lengthen the matrix

# Run function
output <- lengthen_mut_matrix(input)
output_indel <- lengthen_mut_matrix(input_indel)



test_that("Output has correct class", {
  expect_true(inherits(output, "matrix"))
  expect_true(inherits(output_indel, "matrix"))
})

nr_regions <- input %>%
  colnames() %>%
  stringr::str_remove(".*\\.") %>%
  unique() %>%
  length()

nr_regions_indel <- input_indel %>%
  colnames() %>%
  stringr::str_remove(".*\\.") %>%
  unique() %>%
  length()

test_that("Output has correct size", {
  expect_equal(dim(output), c(nrow(input) * nr_regions, ncol(input) / nr_regions))
  expect_equal(dim(output_indel), c(nrow(input_indel) * nr_regions_indel, ncol(input_indel) / nr_regions_indel))
})

expected <- readRDS(system.file("states/mut_mat_longregions.rds",
  package = "MutationalPatterns"
))

expected_indel <- readRDS(system.file("states/blood_indels_longmatrix_split_region.rds",
  package = "MutationalPatterns"
))

test_that("Output transforms correctly", {
  expect_equal(output, expected)
  expect_equal(output_indel, expected_indel)
})
