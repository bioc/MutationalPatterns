context("test-plot_96_profile")


# Load mutation matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Plot the 96-profile of three samples
output <- plot_96_profile(mut_mat[, c(1, 4, 7)])

# Plot a condensed profile
output_condensed <- plot_96_profile(mut_mat[, c(1, 4, 7)], condensed = TRUE)

# Load extracted signatures and plot
nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
  package = "MutationalPatterns"
))
output_signatures <- plot_96_profile(nmf_res$signatures)


# Perform tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_condensed, c("gg")))
  expect_true(inherits(output_signatures, c("gg")))
})
