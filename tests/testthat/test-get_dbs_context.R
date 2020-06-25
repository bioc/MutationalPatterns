context("test-get_dbs_context")


## Get GRangesList with DBS.
dbs_grl <- readRDS(system.file("states/blood_grl_dbs.rds",
  package = "MutationalPatterns"
))

## Set context dbs
output <- get_dbs_context(dbs_grl)

expected <- readRDS(system.file("states/blood_grl_dbs_context.rds",
  package = "MutationalPatterns"
))


test_that("Output has correct class", {
  expect_true(inherits(output, c("GRanges", "CompressedGRangesList")))
})

test_that("Output is equal to expected", {
  expect_equal(output, expected)
})
