context("test-plot_192_profile")

# Load mutation matrix
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
  package = "MutationalPatterns"
))

# Plot profile for some of the samples
output <- plot_192_profile(mut_mat_s[, c(1, 4, 7)])

# You can create a more condensed version of the plot
output_condensed <- plot_192_profile(mut_mat_s[, c(1, 4, 7)], condensed = TRUE)

# Load extracted signatures and plot
nmf_res_strand <- readRDS(system.file("states/nmf_res_strand_data.rds",
  package = "MutationalPatterns"
))
output_signatures <- plot_192_profile(nmf_res_strand$signatures)

# Perform tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_condensed, c("gg")))
  expect_true(inherits(output_signatures, c("gg")))
})
