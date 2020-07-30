context("test-fit_to_signatures")

# Get mut_mat
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Get signatures
signatures <- get_known_signatures()

# Run function
output <- fit_to_signatures(mut_mat, signatures)

# Get expected
expected <- readRDS(system.file("states/snv_refit.rds",
  package = "MutationalPatterns"
))

# Run tests
test_that("Output has correct class", {
  expect_true(inherits(output, "list"))
  expect_true(inherits(output$contribution, "matrix"))
  expect_true(inherits(output$reconstructed, "matrix"))
})

test_that("Output is equal to expected", {
  expect_equal(output, expected)
})

# Get indel mut_mat
indel_counts <- readRDS(system.file("states/blood_indel_counts.rds", package = "MutationalPatterns"))

# Get indel signatures
signatures <- get_known_signatures("indel")

# Get expected
expected <- readRDS(system.file("states/indel_refit.rds",
  package = "MutationalPatterns"
))

# Run tests
test_that("Refitting indels gives expected output.", {
  output <- fit_to_signatures(indel_counts, signatures)
  expect_equal(output, expected)
})

# Get dbs mut_mat
dbs_counts <- readRDS(system.file("states/blood_dbs_counts.rds", package = "MutationalPatterns"))

signatures <- get_known_signatures("dbs")


expected <- readRDS(system.file("states/dbs_refit.rds",
  package = "MutationalPatterns"
))

test_that("Refitting dbss gives expected output.", {
  output <- fit_to_signatures(dbs_counts, signatures)
  expect_equal(output, expected)
})
