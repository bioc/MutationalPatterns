context("test-plot_dbs_contexts")

## Get dbs counts
dbs_counts <- readRDS(system.file("states/blood_dbs_counts.rds",
  package = "MutationalPatterns"
))

## Plot contexts
output <- plot_dbs_contexts(dbs_counts)
output_samey <- plot_dbs_contexts(dbs_counts, same_y = TRUE)
output_condensed <- plot_dbs_contexts(dbs_counts, condensed = TRUE)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_samey, c("gg")))
  expect_true(inherits(output_condensed, c("gg")))
})
