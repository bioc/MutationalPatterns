context("test-plot_original_vs_reconstructed")

# Load mutation matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Load the nmf res
nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
  package = "MutationalPatterns"
))


# Load signature refit.
fit_res <- readRDS(system.file("states/snv_refit.rds",
  package = "MutationalPatterns"
))

# Run function
output <- plot_original_vs_reconstructed(mut_mat, nmf_res$reconstructed)
output_fit <- plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed)
output_intercept <- plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, y_intercept = 0.90)
output_lims <- plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, ylims = c(0, 1))

# Test
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_fit, c("gg")))
  expect_true(inherits(output_intercept, c("gg")))
  expect_true(inherits(output_lims, c("gg")))
})
