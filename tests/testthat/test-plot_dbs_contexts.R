context("test-plot_dbs_contexts")

## Get dbs counts
dbs_counts <- readRDS(system.file("states/blood_dbs_counts.rds",
  package = "MutationalPatterns"
))

## Plot contexts
output <- plot_dbs_contexts(dbs_counts)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
})
