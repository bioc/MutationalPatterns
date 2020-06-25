context("test-mut_type")

# Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))
# Get mut type
output <- mut_type(vcfs[[1]])

# Unit tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("character")))
})

test_that("The 6 base mutation types are returned", {
  base_types <- sort(unique(output))
  expect_equal(base_types, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
})
