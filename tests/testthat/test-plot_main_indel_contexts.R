context("test-plot_main_indel_contexts")


## Get indel counts
indel_counts <- readRDS(system.file("states/blood_indel_counts.rds",
  package = "MutationalPatterns"
))

## Plot contexts
output <- plot_main_indel_contexts(indel_counts)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
})
