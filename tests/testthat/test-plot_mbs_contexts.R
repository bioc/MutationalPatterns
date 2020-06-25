context("test-plot_mbs_contexts")

## Get mbs counts
mbs_counts <- readRDS(system.file("states/blood_mbs_counts.rds",
  package = "MutationalPatterns"
))

## Plot contexts
output <- plot_mbs_contexts(mbs_counts)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
})
