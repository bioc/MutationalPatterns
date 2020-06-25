context("test-plot_compare_mbs")

# Get the mbs counts
mbs_counts <- readRDS(system.file("states/blood_mbs_counts.rds",
  package = "MutationalPatterns"
))


# Run default function
output <- plot_compare_mbs(
  mbs_counts[, 1],
  mbs_counts[, 2]
)

# Change the names of the profiles
output_name <- plot_compare_mbs(mbs_counts[, 1],
  mbs_counts[, 2],
  profile_names = c("Original", "Reconstructed")
)

# Change the y_limits
output_yaxis <- plot_compare_mbs(mbs_counts[, 1],
  mbs_counts[, 2],
  profile_ymax = 0.9,
  diff_ylim = c(-0.8, 0.8)
)

# Perform tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_name, c("gg")))
  expect_true(inherits(output_yaxis, c("gg")))
})
