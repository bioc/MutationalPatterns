context("test-plot_compare_indels")

# Get indel counts
indel_counts <- readRDS(system.file("states/blood_indel_counts.rds",
  package = "MutationalPatterns"
))

# Get indel refit
fit_res <- readRDS(system.file("states/indel_refit.rds",
  package = "MutationalPatterns"
))

# Run default function
output <- plot_compare_indels(indel_counts[, 1], fit_res$reconstructed[, 1])

# Test you can change the name
output_name <- plot_compare_indels(indel_counts[, 1],
  fit_res$reconstructed[, 2],
  profile_names = c("Original", "Reconstructed")
)

## You can also change the y limits.
## This can be done separately for the profiles and the different facets.
output_yaxis <- plot_compare_indels(indel_counts[, 1],
  fit_res$reconstructed[, 2],
  profile_ymax = 0.3,
  diff_ylim = c(-0.03, 0.03)
)

# Perform tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_name, c("gg")))
  expect_true(inherits(output_yaxis, c("gg")))
})
