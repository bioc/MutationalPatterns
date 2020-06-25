context("test-count_dbs_contexts")

## Get a GRangesList object with dbs contexts.
grl_dbs_context <- readRDS(system.file("states/blood_grl_dbs_context.rds",
  package = "MutationalPatterns"
))

output <- count_dbs_contexts(grl_dbs_context)
expected <- readRDS(system.file("states/blood_dbs_counts.rds",
  package = "MutationalPatterns"
))

test_that("Output has correct class", {
  expect_true(inherits(output, c("matrix")))
})

test_that("Output is identical to expected", {
  expect_identical(output, expected)
})
