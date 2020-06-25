context("test-plot_compare_dbs")

# Get dbs counts
dbs_counts <- readRDS(system.file("states/blood_dbs_counts.rds",
  package = "MutationalPatterns"
))

# Get dbs refit
fit_res <- readRDS(system.file("states/dbs_refit.rds",
  package = "MutationalPatterns"
))

# Run default function
output <- plot_compare_dbs(dbs_counts[, 1], fit_res$reconstructed[, 1])

# Test you can change the name
output_name <- plot_compare_dbs(dbs_counts[, 1],
  fit_res$reconstructed[, 2],
  profile_names = c("Original", "Reconstructed")
)

## You can also change the y limits.
## This can be done separately for the profiles and the different facets.
output_yaxis <- plot_compare_dbs(dbs_counts[, 1],
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
