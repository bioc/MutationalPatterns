context("test-plot_indel_contexts")

## Get indel counts
indel_counts <- readRDS(system.file("states/blood_indel_counts.rds",
  package = "MutationalPatterns"
))

## Plot contexts
output <- plot_indel_contexts(indel_counts)
output_same_y <- plot_indel_contexts(indel_counts, same_y = TRUE)
output_extra_labels <- plot_indel_contexts(indel_counts, extra_labels = TRUE)
output_condensed <- plot_indel_contexts(indel_counts, condensed = TRUE)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_same_y, c("gg")))
  expect_true(inherits(output_extra_labels, c("gg")))
  expect_true(inherits(output_condensed, c("gg")))
})
