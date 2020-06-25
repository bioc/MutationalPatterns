context("test-count_mbs_contexts")

## Get a GRangesList object with mbs contexts.
grl_mbs <- readRDS(system.file("states/blood_grl_mbs.rds",
  package = "MutationalPatterns"
))

output <- count_mbs_contexts(grl_mbs)
expected <- readRDS(system.file("states/blood_mbs_counts.rds",
  package = "MutationalPatterns"
))

test_that("Output has correct class", {
  expect_true(inherits(output, c("matrix")))
})

test_that("Output is identical to expected", {
  expect_identical(output, expected)
})
