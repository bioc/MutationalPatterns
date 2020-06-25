context("test-count_indel_contexts")

## Get a GRangesList object with indel contexts.
grl_indel_context <- readRDS(system.file("states/blood_grl_indel_context.rds",
  package = "MutationalPatterns"
))

output <- count_indel_contexts(grl_indel_context)
expected <- readRDS(system.file("states/blood_indel_counts.rds",
  package = "MutationalPatterns"
))

test_that("Output has correct class", {
  expect_true(inherits(output, c("matrix")))
})

test_that("Output is identical to expected", {
  expect_identical(output, expected)
})
